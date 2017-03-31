#include "type.h"

module ot_expr
  use ot_common
  use ot_tensor
  use dispmodule
  use ot_ref
  use ot_vector
  use ot_geom
  use ot_node

!#define ASSERT(x,msg) call assert(x, __FILE__, __LINE__, msg)
#:set DEBUG = 1
#define DEBUG
#:include "type_def.fypp"
  
  public :: operator(+),operator(-)
  public :: operator(*),operator(/)
  public :: assignment(=)

  interface assignment (=)
     module procedure tensor_assign_tensor
     module procedure node_assign_tensor
     module procedure node_assign_node
  end interface assignment (=)


  interface slice
#:set slice_idx_type=['int', 'arr','range', 'char']
#:for data in ['node', 'tensor']
#:for ta in slice_idx_type
     module procedure slice_${data}$_${ta}$     
#:for tb in slice_idx_type
     module procedure slice_${data}$_${ta}$_${tb}$     
#:for tc in slice_idx_type
     module procedure slice_${data}$_${ta}$_${tb}$_${tc}$
#:endfor
#:endfor
#:endfor
#:endfor
  end interface slice

  interface set
#:for src_type in ['tensor', 'node', 'int', 'real', 'real8']
#:for dst_type in ['tensor', 'node']
     module procedure set_${dst_type}$_${src_type}$
#:endfor
#:endfor
  end interface set
 
  interface operator (**)
#:for type1 in ['node', 'tensor']     
#:for type2 in ['integer', 'real', 'real8']
     module procedure ${type1}$_pow_${type2}$
#:endfor
#:endfor
  end interface operator (**)
  
  !the following code using preprossor to create interfaces
#:for op in [['plus','+'], ['minus','-'], ['mult','*'], ['divd','/'], &
  ['gt','>'], ['ge', '>='], ['lt', '<'],['le', '<='],['eq','=='],['ne','/=']]
  interface operator (${op[1]}$)
#:for type1 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:for type2 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:if any(type1 in s for s in ['real(8)', 'real', 'integer']) & 
     and any(type2 in s for s in ['real(8)', 'real', 'integer'])
#:else
#:set name1 = re.sub('\(8\)', '8', type1)
#:set name2 = re.sub('\(8\)', '8', type2)
     module procedure ${name1}$_${op[0]}$_${name2}$
#:endif
#:endfor
#:endfor
  end interface operator (${op[1]}$)
  
#:endfor

!interface for exp, log, pow, ...
#:for e in L
#:if e[1] >= 111
#:if e[3] == '+' or e[3] == '-'
  interface operator(${e[3]}$)
#:else
  interface ${e[3]}$
#:endif
#:for type1 in ['tensor', 'node']
     module procedure ${e[2]}$_${type1}$
#:endfor
  end interface

#:endif
#:endfor

  interface disp
     module procedure disp_node
  end interface disp

contains
  
  subroutine init_expr(ierr)
    use ot_tensor
    implicit none
    integer, intent(out) :: ierr
  end subroutine

  subroutine disp_node(A, msg)
    implicit none
    type(node) :: A
    type(tensor) :: T
    character(len=*),intent(in),optional :: msg
    integer :: ierr
    
    T = A !this will evaluate the node
    if(present(msg)) then
       call disp(T, msg)
    else
       call disp(T)
    end if
    call destroy(T, ierr)
  end subroutine

  subroutine node_assign_node(A, B)
    implicit none
    type(node), intent(inout), target :: A
    type(node), intent(in),    target :: B

#:if DEBUG > 0
    print*, "call node_assign_node"
#:endif
    
    if(associated(B%data)) &
         call assign_ptr(A%data, B%data)
    
    A%node_type = B%node_type
    A%m_shape = B%m_shape
    A%alpha = B%alpha
    A%beta = B%beta
    A%scalar = B%scalar
    call copy_node_ptr(A%operands, B%operands)
    
  end subroutine
  
  subroutine node_assign_tensor(A, B)
    implicit none

    type(tensor), intent(out) :: A
    
    !the type of B must be "intent(in)"
    type(node), intent(in)   :: B 
    type(node) :: C
    type(node_ptr), allocatable :: res(:)
    integer ierr
    
    call node_new(C, B)

    ! call disp(B%ref%ref_box, 'B%ref_box = ')
    ! call disp(C%ref%ref_box, 'C%ref_box = ')       

    !call node_optimize(C)

    !call disp_tree(C)
    
    !call disp_info(C, 'C=')
    
    ! write(*, "(A, Z16.16)"), "op1=", loc(C%operands(1)%ptr)
    ! write(*, "(A, Z16.16)"), "op2=", loc(C%operands(2)%ptr)
    
    call write_graph(C, file='C.dot')
    
    call eval(A, C, ierr, .true.)

    ! call node_destroy(C, ierr)
    
    !call tensor_copy(A, B%data)
    
    !if(B%node_type /= type_data) call eval(A)
    ! if(B%node_type == type_data) then
    !    call tensor_copy(A, B)
    ! else
    !    call tensor_copy(C, B)
    !    call eval(C)
    !    call tensor_copy(A, C)
    ! endif
    
    if(allocated(res)) deallocate(res)    
  end subroutine

  !*********************************
  !> evaluate expression nodes
  !*********************************
  recursive subroutine eval(res, A, ierr, force_eval_arg)
    implicit none
    type(tensor), intent(inout) :: res
    type(node), intent(in) :: A
    integer :: i, num_oprands
    integer, intent(out) :: ierr
    real(8), allocatable :: ops_alpha(:), ops_beta(:)
    real(8) :: alpha, beta
    type(tensor_ptr), allocatable :: ops_tensor(:)
    logical, optional, intent(in) :: force_eval_arg
    logical :: force_eval
    type(node), pointer :: node_ptr
    
    ! interface
    !    subroutine invoke(p, res, operands, alpha, beta, args, n) &
    !         bind(C, name="invoke")
    !      implicit none
    !      integer :: p
    !      C_POINTER :: res
    !      C_POINTER :: operands
    !      C_POINTER :: alpha
    !      C_POINTER :: beta
    !      C_POINTER :: args
    !      integer :: n
    !    end subroutine
    ! end interface

    if(present(force_eval_arg)) then
       force_eval = force_eval_arg
    else
       force_eval = .true.
    endif

    !call disp_info(A, 'A = ')

    !process the reference node
    if(is_ref(A)) then
       if(allocated (A%operands)) then
          call assert(size(A%operands) == 1, &
               __FILE__, __LINE__, &
               "the number of operands &
               for ref node can not be larger than 1")
          call eval(A%data, A%operands(1)%ptr, ierr, .false.)
       end if

       call slice_tensors(res, A%data, A%ref)
       return
    end if

    !if this node is the root node, force to evaluate
    if(is_data(A) .and. .not. force_eval) return
    
    !if(allocated(A%operands)) then
    if(.not. is_data(A)) then       
       num_oprands = size(A%operands)

       allocate(ops_tensor(num_oprands), &
            ops_alpha(num_oprands), ops_beta(num_oprands))
       
       !recusively evaluate sub-nodes
       do i = 1, num_oprands

          ! call disp_info(A, 'A = ')
          ! call disp_info(A%operands(i)%ptr, 'A1 = ')
          
          node_ptr => A%operands(i)%ptr

          if(.not. associated(node_ptr%data)) then
             allocate(node_ptr%data)
          endif
          
          call eval(node_ptr%data, node_ptr, ierr, .false.)

          ! print*, "============================"
          ! call disp_info(node_ptr, 'NODE')
          ! call disp(node_ptr%data, 'TENSOR = ')

          call assign_ptr(ops_tensor(i)%ptr, node_ptr%data)
          ops_alpha(i) = node_ptr%alpha
          ops_beta(i)  = node_ptr%beta
       end do
    else
       allocate(ops_tensor(1), ops_alpha(1), ops_beta(1))
       call assign_ptr(ops_tensor(1)%ptr, A%data)
       ops_alpha(1) = A%alpha
       ops_beta(1)  = A%beta
    endif

    call tensor_duplicate(res, ops_tensor(1)%ptr)
    
    !call invoke(A%node_type, loc(A%data), loc(ops_tensor), &
    !     loc(alpha), loc(beta), loc(A%args), num_oprands)
    !write(*, "(A, Z16.6)") "A2=", loc(A%data)

    !if(need_eval(A) .or. force_eval) then
       select case(A%node_type)
#:for e in L
#:if e[4]=='A' or e[4]=='B' or e[4]=='C'
       case (${e[1]}$)
#ifdef DEBUG
          print*, "calling function ${e[2]}$ ..."
#endif
          call ${e[2]}$_tensors(res, &
               ops_tensor, ops_alpha, ops_beta, A%args)
#:endif
#:endfor
       case default
          print*, "calling function eval_tensors ..."          
          call eval_tensors(res, &
               ops_tensor, ops_alpha, ops_beta, A%args)
       end select

       ! if(.not. associated(A%data)) then
       !    print*, "execute function ", op_names(A%node_type), " failed."
       !    call abort()
       ! end if
    !endif
      
    print*, "evalation of ", trim(op_names(A%node_type)), " is complete."
  end subroutine
  
  !> this function will be called in the following two situations:
  !> 1) A = B, where A and B are both tensors
  !> 2) func(A), where A is a tensor and passed to a function as an argument.
  !> in both case, we do not need a deep copy
  subroutine tensor_assign_tensor(A, B)
    implicit none
    type(tensor), intent(out) :: A
    type(tensor), intent(in)  :: B

    print*, "--> called tensor_assign_tensor"
    
    call tensor_deep_copy(A, B)
  end subroutine
  

  recursive subroutine node_optimize(A) 
    implicit none
    type(node), intent(inout), target :: A
    type(node_ptr),  allocatable :: subs(:), tmp(:)
    integer i, cnt, num_operands, j
    type(node), pointer :: child
    
    if(is_data(A)) return

    num_operands = size(A%operands)

    call disp_info(A)

    i = 1
    do while(i <= num_operands)
       child => A%operands(i)%ptr
       if(A%node_type == child%node_type &
            .and. is_arithmetic(A)) then

          call remove(A%operands, i)
          call push_back(A%operands, child%operands)

          num_operands = size(A%operands)
       else
          i = i + 1
       end if
    enddo

    !call disp_info(A)
    
    do i = 1, size(A%operands) 
       call node_optimize(A%operands(i)%ptr)
    enddo

  end subroutine
  
  recursive subroutine write_graph(A, is_root, file)
    implicit none
    type(node), intent(in) :: A
    integer :: i
    character(len=100) :: label1, label2
    character(len=*),optional :: file
    integer id
    logical,intent(in), optional :: is_root
    logical :: am_i_root
    integer, parameter :: out_unit = 20
    character(len=10) :: op_name
    character(len=*), parameter :: &
         format_label = "(I0.1,A,A,A,I0.1,A,F4.1,A,F4.1,A)"
    character(len=*), parameter :: &
         format_label1 = "(I0.1,A,A,A)"
    character(len=*), parameter :: &
         format_map = "(I0.1,A,I0.1,A)"

    !only execute this function on first process
    if(get_rank() > 0) return
    
    if(present(is_root)) then
       am_i_root = is_root
    else
       am_i_root = .true.
    endif

    if(am_i_root) then
       call reset_global_id()
       if(present(file)) then
          open(unit=out_unit,file=file, &
               action="write", status="replace")
       else
          open(unit=out_unit,file="graph.dot", &
               action="write", status="replace")
       endif
       write(out_unit, *) "digraph G {"
    end if

    write(op_name, *) trim(op_names(A%node_type))

    ! id = get_global_id()
    ! A%id = id
    
#:set SHOW_ALPHA_BETA = 1
    
#:if SHOW_ALPHA_BETA > 0
    write(out_unit, format_label) A%id, &
         '[label="[',trim(adjustl(op_name)),']\n id=', A%id, '\n(',&
         A%alpha, ',',A%beta, ')"];'
#:else
    write(out_unit, format_label1) A%id, &
         '[label="',trim((adjustl(op_name)), '"];'
#:endif
    
    !print*, label1
    if(is_data(A)) then
       if(am_i_root) then
          write(out_unit, *) "}"
          close (out_unit)
       end if
       return
    end if
    
    if(allocated(A%operands)) then
       do i = 1, size(A%operands)
          call write_graph(A%operands(i)%ptr,  .false.)
          write(out_unit, format_map), &
               A%id, "->", A%operands(i)%ptr%id, ";"
       end do
    endif
    
    if(am_i_root) then
       write(out_unit, *) "}"
       close (out_unit)
    end if
    
  end subroutine

  !the following code using preprossor to create subroutines  
#:for type1 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:for type2 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:for op in ['plus', 'minus', 'mult', 'divd', 'gt', 'ge', 'lt', 'le', 'eq', 'ne']
#:if any(type1 in s for s in ['real(8)', 'real', 'integer']) & 
     and any(type2 in s for s in ['real(8)', 'real', 'integer'])  
#:else  
#:set name1 = re.sub('\(8\)', '8', type1)
#:set name2 = re.sub('\(8\)', '8', type2)
  function ${name1}$_${op}$_${name2}$ (A, B) result(res)
    implicit none       
    type(${type1}$),intent(in),target  :: A
    type(${type2}$),intent(in),target  :: B
    type(node), target :: res
    type(node), pointer :: C, D

    allocate(C, D)

    call node_new(C, A)
    call node_new(D, B)
    call node_new(res, type_${op}$, C, D)
  end function
  
#:endif
#:endfor
#:endfor
#:endfor

#:for e in L
#:if e[1] >= 111
  #:for type in ['tensor', 'node']
  !> function ${e[2]}$(${type}$), return a tree node 
  function ${e[2]}$_${type}$ (o) result(res)
    implicit none
    type(${type}$), intent(in) :: o
    type(node) :: res

    call node_new(res, ${e[0]}$, o)
    
  end function

#:endfor
#:endif
#:endfor

#:for type1 in ['node', 'tensor']
#:for type2 in ['integer', 'real', 'real(8)']
#:set name2 = re.sub('\(8\)', '8', type2)
  !> power function ${type1}$**${type2}$
  function ${type1}$_pow_${name2}$(o, n) result(res)
    implicit none
    type(${type1}$), intent(in) :: o
    type(${type2}$), intent(in) :: n
    type(node) :: res
    type(node), pointer :: new_o

    allocate(new_o)
    call node_new(new_o, o)
    call node_new(res, type_pow, new_o)
    res%args(1) = n
    
  end function
  
#:endfor
#:endfor


  !*****************************
  !> slice functions
  !*****************************
#:set itype=[['int', 'integer'], &
       ['range', 'type(range)'], &
       ['arr', 'integer,dimension(:)'], &
       ['char', 'character(len=1)']]
       
#:for data in ['node', 'tensor']  
#:for ta in itype
   function slice_${data}$_${ta[0]}$ &
        (obj, a) result(res)
    implicit none
    type(${data}$),target :: obj    
    ${ta[1]}$ :: a
    type(node), pointer :: res, tmp
    integer :: dim
    integer :: xs, xe, ys, ye, zs, ze
    integer :: dim_x, dim_y, dim_z

    res => slice_${data}$_${ta[0]}$_int_int(obj, a, 1, 1)
  end function

#:for tb in itype
  function slice_${data}$_${ta[0]}$_${tb[0]}$ &
       (obj, a, b) result(res)
    implicit none
    type(${data}$),target :: obj    
    ${ta[1]}$ :: a
    ${tb[1]}$ :: b    
    type(node), pointer :: res, tmp
    integer :: dim
    integer :: xs, xe, ys, ye, zs, ze
    integer :: dim_x, dim_y, dim_z

    res => slice_${data}$_${ta[0]}$_${tb[0]}$_int(obj, a, b, 1)
  end function
       
#:for tc in itype
  function slice_${data}$_${ta[0]}$_${tb[0]}$_${tc[0]}$ &
       (obj, a, b, c) result(res)
    implicit none
    type(${data}$),target :: obj    
    ${ta[1]}$ :: a
    ${tb[1]}$ :: b    
    ${tc[1]}$ :: c
    type(node), pointer :: res, tmp
    integer :: dim
    integer :: xs, xe, ys, ye, zs, ze
    integer :: dim_x, dim_y, dim_z
    integer :: src_shape(3)

    allocate(res)

    call assert(is_valid(obj), &
         __FILE__, __LINE__, &
         'the input tensor/node must has a valid shape')
    
    src_shape = shape(obj)

#:if ta[0] == 'int'
    dim_x = 1
    res%ref%range_x = r(a-1, a-1)
    res%ref%ref_index_type_x = 0
    call assert(a>=1 .and. a<=src_shape(1), &
         __FILE__, __LINE__, &
         "x index exceeds array dimensions")
#:elif ta[0] == 'char'
    dim_x = src_shape(1)
    res%ref%range_x = r(0, dim_x-1)
    res%ref%ref_index_type_x = 0
#:elif ta[0] == 'range'
    if(a == range_all) then
       dim_x = src_shape(1)
       res%ref%range_x = r(0, dim_x - 1)
    else 
       dim_x = a%upper - a%lower + 1
       res%ref%range_x = a - 1
    endif
    res%ref%ref_index_type_x = 0
    call assert(a%lower>=1 .and. a%upper<=src_shape(1), &
         __FILE__, __LINE__, &
         "x index exceeds array dimensions")    
#:else
    allocate(res%ref%iarr_x(size(a)))
    res%ref%iarr_x = a - 1
    dim_x = size(a)
    res%ref%ref_index_type_x = 1    
#:endif

#:if tb[0] == 'int'
    dim_y = 1
    res%ref%range_y = r(b-1,b-1)
    res%ref%ref_index_type_y = 0
    call assert(b>=1 .and. b<=src_shape(2), &
         __FILE__, __LINE__, &
         "y index exceeds array dimensions")
#:elif tb[0] == 'char'
    dim_y = src_shape(2)
    res%ref%range_y = r(0, dim_y-1)
    res%ref%ref_index_type_y = 0    
#:elif tb[0] == 'range'
    if(b == range_all) then
       dim_y = src_shape(2)
       res%ref%range_y = r(0, dim_y-1)
    else 
       dim_y = b%upper - b%lower + 1
       res%ref%range_y = b-1
    endif
    
    res%ref%ref_index_type_y = 0
    call assert(b%lower>=1 .and. b%upper<=src_shape(2), &
         __FILE__, __LINE__, &
         "y index exceeds array dimensions")        
#:else
    allocate(res%ref%iarr_y(size(b)))
    res%ref%iarr_y = b-1
    dim_y = size(b)
    res%ref%ref_index_type_y = 1    
#:endif

#:if tc[0] == 'int'
    dim_z = 1
    res%ref%range_z = r(c-1, c-1)
    res%ref%ref_index_type_z = 0
    call assert(c>=1 .and. c<=src_shape(3), &
         __FILE__, __LINE__, &
         "z index exceeds array dimensions")
#:elif tc[0] == 'char'
    dim_z = src_shape(3)
    res%ref%range_z = r(0, dim_z-1)
    res%ref%ref_index_type_z = 0
#:elif tc[0] == 'range'
    if(c == range_all) then
       dim_z = src_shape(3)
       res%ref%range_z = r(0, dim_z - 1)
    else 
       dim_z = c%upper - c%lower + 1
       res%ref%range_z = c-1
    endif
    res%ref%ref_index_type_z = 0
    call assert(c%lower>=1 .and. c%upper<=src_shape(3), &
         __FILE__, __LINE__, &
         "z index exceeds array dimensions")            
#:else
    allocate(res%ref%iarr_z(size(c)))
    res%ref%iarr_z = c-1
    dim_z = size(c)
    res%ref%ref_index_type_z = 1    
#:endif
    
#:if data == 'tensor'
    !call node_new(res, obj) !create a data node
    call assign_ptr(res%data, obj)
    !call node_add_operand(res, tmp)
#:else
    call assign_ptr(tmp, obj) 
    call node_add_operand(res, tmp)
#:endif

    res%m_shape = (/dim_x, dim_y, dim_z/)
    res%node_type = type_ref
    
  end function
  
#:endfor
#:endfor
#:endfor
#:endfor  

#:for src_type in [['tensor', 'tensor', 'arr'], &
  ['node', 'node', 'arr'], ['int', 'integer', 'scalar'], &
       ['real', 'real', 'scalar'], ['real8', 'real(8)', 'scalar']]
#:for dst_type in [['tensor', 'tensor'], ['node', 'node']]
  
  subroutine set_${dst_type[0]}$_${src_type[0]}$(dst, src)
    type(${dst_type[1]}$), intent(in) :: dst
    type(${src_type[1]}$) :: src
    type(ref_info) :: set_ref

#:if src_type[2] == 'scalar' and dst_type[0]=='tensor' 
    call range_to_ref (set_ref, range_all, range_all, range_all)
    call data_set_scalar (dst%data, real(src, 8), set_ref)
#:endif
    
#:if src_type[2] == 'scalar' and dst_type[0] == 'ref' 
    call data_set_scalar (dst%data%data, real(src, 8),  dst%ref)
#:endif

#:if src_type[2] == 'tensor' and dst_type[0] == 'tensor' 
#:endif

#:if src_type[2] == 'tensor' and dst_type[0] == 'node'
    call assert(is_ref(dst), __FILE__, __LINE__,&
         "the target must be a reference.")
    call set_node_node(dst, slice(src, range_all, range_all, range_all))
#:endif

#:if src_type[1] == 'node' and dst_type[1] == 'node'
    call assert(associated(dst%data), __FILE__, __LINE__,&
         "as lvalue, node%data must be associated.")

    call assert(is_ref(dst), __FILE__, __LINE__,&
         "the target must be a reference.")
    
    if(allocated(src%operands) .and. &  !for expression node    
         size(src%operands)>1) then

       !if the source node is an expression node,
       !evaluate it firt
       call eval(src%data, src, ierr)
       
       call data_set_ref_ref(dst%data%data, &
            src%data%data, dst%ref, src%ref)
       
    else if(src%node_type == type_data) then !for data node    
       !if the source node is an data node,
       !we have to create reference info using its shape
       call data_set_ref_ref(dst%data%data, &
            src%data%data, shape_to_ref(shape(src)), &
            src%ref)
    else
       call data_set_ref_ref(dst%data%data, &
            src%data%data, dst%ref, src%ref)
    end if
#:endif
    
  end subroutine
#:endfor
#:endfor

end module

