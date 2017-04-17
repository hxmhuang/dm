#include "type.h"

#define MAX_TREESTR_LENGTH 100000
#define SKIP_EW_NODE 1
#define ALL_NODE 1

module ot_expr
  use ot_common
  use ot_array
  use dispmodule
  use ot_ref
  use ot_geom
  use ot_node
  use ot_kernels
  use ot_buffer
  
!#define ASSERT(x,msg) call assert(x, __FILE__, __LINE__, msg)
#:set DEBUG = 1
#define DEBUG
#:include "type_def.fypp"
  
  public :: operator(+),operator(-)
  public :: operator(*),operator(/)
  public :: assignment(=)

  interface assignment (=)
     module procedure array_assign_array
     module procedure node_assign_array
     module procedure node_assign_node
#:for type in ['i4', 'r4', 'r8']
     module procedure array_assign_${type}$
#:endfor
  end interface assignment (=)

  interface disp
     module procedure disp_node
  end interface disp

  interface slice
#:set slice_idx_type=['int', 'arr','range', 'char']
#:for data in ['node', 'array']
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
#:for src_type in ['array', 'node', 'int', 'real', 'real8']
#:for dst_type in ['array', 'ref']
     module procedure set_${dst_type}$_from_${src_type}$
#:endfor
#:endfor
  end interface set
 
  !the following code using preprossor to create interfaces
#:for op in [['plus','+'], ['minus','-'], ['mult','*'], ['divd','/'], &
  ['gt','>'], ['ge', '>='], ['lt', '<'],['le', '<='],['eq','=='],['ne','/=']]
  interface operator (${op[1]}$)
#:for type1 in ['real(8)', 'real', 'integer', 'array', 'node']
#:for type2 in ['real(8)', 'real', 'integer', 'array', 'node']
#:if any(type1 in s for s in ['real(8)', 'real', 'integer']) & 
     and any(type2 in s for s in ['real(8)', 'real', 'integer'])
#:else
#:set name1 = re.sub('\(8\)', '8', type1)
#:set name2 = re.sub('\(8\)', '8', type2)
     module procedure expr_${name1}$_${op[0]}$_${name2}$
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
#:for type1 in ['array', 'node']
     module procedure expr_${e[2]}$_${type1}$
#:endfor
  end interface

#:endif
#:endfor

  interface operator (**)
#:for type1 in ['node', 'array']     
#:for type2 in ['integer', 'real', 'real8']
     module procedure expr_${type1}$_pow_${type2}$
#:endfor
#:endfor
  end interface operator (**)

  logical :: is_opt_mode = .false. ! if mode=1, the expression just generate kernels with evaluation
  character(len=*), parameter :: kernel_file = "kernels.fypp"
contains

  subroutine init_expr(ierr)
    use ot_array
    implicit none
    integer, intent(out) :: ierr
    integer :: out_unit = 10
    
    if(ot_has_option('--opt_mode')) then
       is_opt_mode = .true.
       !call ot_option_int('--opt_mode', opt_mode, ierr)
       if(get_rank(MPI_COMM_WORLD) == 0) then
          print*, "===================================================="
          print*, "INFOMATION : Program is runing in optimization mode, &
               & only optimization info generated. &
               & Evaluation will NOT be performed!"
          print*, "===================================================="

          !clean up the kernel file
          open(unit=out_unit,file=kernel_file, &
               action='write', status='replace')
          close(out_unit)
          
       end if
    endif

  end subroutine

  subroutine disp_node(A, msg)
    implicit none
    type(node), intent(inout), pointer :: A
    character(len=*),intent(in),optional :: msg
    type(array) :: T    
    integer :: ierr

    T = A !this will evaluate the node
    
    if(present(msg)) then
       call disp(T, msg)
    else
       call disp(T)
    end if

    call destroy(T, ierr)

    if(is_rvalue(A)) then
       call destroy(A, ierr)
    end if 
  end subroutine
  
  subroutine node_assign_node(A, B)
    implicit none
    type(node), intent(inout) :: A
    type(node), intent(in) :: B

! #:if DEBUG > 0
!     print*, "call node_assign_node"
! #:endif
    
    if(associated(B%data)) &
         call assign_ptr(A%data, B%data)
    
    A%node_type = B%node_type
    A%m_shape = B%m_shape
    A%alpha = B%alpha
    A%beta = B%beta
    A%scalar = B%scalar
    call copy_node_ptr(A%operands, B%operands)
    
  end subroutine
  
  subroutine node_assign_array(A, B)
    implicit none
    type(array), intent(inout) :: A
    type(node), intent(in) :: B
    type(node), pointer :: src
    type(array), pointer :: dst
    type(array) :: res
    integer ierr, cnt
    character(len=MAX_TREESTR_LENGTH) :: str
    integer, parameter :: out_unit = 10
    
    ! print*, "node assigning to array ..."

    call assign_ptr(dst, A)
    call assign_ptr(src, B)
    
    ! call disp_info(B, 'B = ')
    ! call destroy(p, ierr)
    ! return
    
    !call node_new(C, B)

    ! call disp(B%ref%ref_box, 'B%ref_box = ')
    ! call disp(C%ref%ref_box, 'C%ref_box = ')       

    !call node_optimize(C)
    !call disp_info(C, 'C=')
    
    ! write(*, "(A, Z16.16)"), "op1=", loc(C%operands(1)%ptr)
    ! write(*, "(A, Z16.16)"), "op2=", loc(C%operands(2)%ptr)
    
    !call write_graph(C, file='C.dot')
    ! print*, '-------------------------'
    ! call disp_tree(src)
    
    if(is_opt_mode) then
       call gen_kernels(src)
    else
       ! print*, "cnt2 = ", expr_sub_cnt(src, 1)    
       ! print*, "cnt1 = ", expr_sub_cnt(src)

       !print*, "src%node_type = ", src%node_type
       !call write_graph(src, file='src.dot')

       ! print*, '-------------------------'           
       ! call disp_tree(src)           
       !call disp_info(src, 'src = ')
       
       call gen_hash(src, .true.)

       ! call disp_info(src)
       ! call disp(src, 'src = ')
       
       call eval(dst, src, ierr, .true.)
       
       ! call disp(dst, 'dst = ') 

       ! dst%m_shape = (/1,2,3/)
       !A%m_shape = (/1,2,3/)       
       ! call disp(A, 'D1 = ')
       
       call set_lvalue(dst)
    endif
    
    call release_ptr(src)
    call release_ptr(dst)
    
    !call disp_tree(B)
    
    ! call node_destroy(C, ierr)
    ! call array_copy(A, B%data)
    
    !if(B%node_type /= type_data) call eval(A)
    ! if(B%node_type == type_data) then
    !    call array_copy(A, B)
    ! else
    !    call array_copy(C, B)
    !    call eval(C)
    !    call array_copy(A, C)
    ! endif
    
    !if(allocated(res)) deallocate(res)    
  end subroutine

  recursive subroutine eval_non_ew(o, ops, alphas, betas, args)
    implicit none
    type(node), intent(inout) :: o
    integer :: i
    integer :: ierr
    
    type(BUFFER_CPTR), intent(inout) :: ops(:)
    type(buffer_r8), intent(inout) :: alphas
    type(buffer_r8), intent(inout) :: betas
    type(buffer_r8), intent(inout) :: args
    
    if(.not. is_ew(o)) then
       call eval(o%data, o, ierr, .true.)
       call push_back(alphas, o%alpha)
       call push_back(betas,  o%beta)
       call push_back(args,   o%args)
    end if
    
    do i = 1, size(o%operands)
       call eval_non_ew(o%operands(i)%ptr, ops, alphas, betas, args)
    enddo
  end subroutine
  
  !*********************************
  !> evaluate expression nodes
  !*********************************
  recursive subroutine eval(res, A, ierr, force_eval_arg)
    implicit none
#include "petsc.h"
    type(array), intent(inout), pointer :: res
    type(node), intent(inout) :: A
    integer :: i, num_oprands
    integer, intent(out) :: ierr
    real(8), allocatable :: ops_alpha(:), ops_beta(:)
    type(array_ptr), allocatable :: ops_array(:)
    logical, optional, intent(in) :: force_eval_arg
    logical :: force_eval
    type(node), pointer :: node_ptr
    type(array), pointer :: tmp_res
    type(buffer_list_r8) :: buf_list
    type(buffer_r8) :: buf
    C_POINTER :: kernel_addr
    type(buffer_r8), pointer :: buf_ptr(:)

    if(present(force_eval_arg)) then
       force_eval = force_eval_arg
    else
       force_eval = .true.
    endif

    !call disp_info(A, 'A = ')
    !> this node correspond to a kernel
    if(A%hash /= 0) then
       !call eval_non_ew(A)
    end if
    
    !process the reference node
    if(is_ref(A)) then
       if(allocated (A%operands)) then
          call assert(size(A%operands) == 1, &
               __FILE__, __LINE__, &
               "the number of operands &
               for ref node can not be larger than 1")
          call eval(A%data, A%operands(1)%ptr, ierr, .false.)
       end if

       call slice_array(res, A%data, A%ref)
       return
    end if

    !if this node is the root node, force to evaluate
    if(is_data(A) .and. .not. force_eval) return

    !call disp_info(A, 'A1 = ')
    
    !if(allocated(A%operands)) then
    if(.not. is_data(A)) then       
       num_oprands = size(A%operands)

       allocate(ops_array(num_oprands), &
            ops_alpha(num_oprands), ops_beta(num_oprands))
       
       !recusively evaluate sub-nodes
       do i = 1, num_oprands

          ! call disp_info(A, 'A = ')
          ! call disp_info(A%operands(i)%ptr, 'A1 = ')
          
          node_ptr => A%operands(i)%ptr

          call eval(node_ptr%data, node_ptr, ierr, .false.)

          ! print*, "============================"
          ! call disp_info(node_ptr, 'NODE')
          ! call disp(node_ptr%data, 'TENSOR = ')

          !call assign_ptr(ops_array(i)%ptr, node_ptr%data)
          ops_array(i)%ptr => node_ptr%data
          ops_alpha(i) = node_ptr%alpha
          ops_beta(i)  = node_ptr%beta
       end do
    else
       allocate(ops_array(1), ops_alpha(1), ops_beta(1))
       !call assign_ptr(ops_array(1)%ptr, A%data)
       ops_array(1)%ptr => A%data
       ops_alpha(1) = A%alpha
       ops_beta(1)  = A%beta
    endif

    !call disp_info(A, 'A2 = ')
    
    if(.not. associated(res)) then
       do i = 1, size(ops_array) 
          if(is_rvalue(ops_array(i)%ptr)) then
             call assign_ptr(res, ops_array(i)%ptr)
             exit
          endif
       end do
    endif
    
    if(.not. associated(res)) then       
       allocate(tmp_res)
       !allocate space for result array object
       call array_duplicate(tmp_res, ops_array(1)%ptr)
       call assign_ptr(res, tmp_res)
       !the evaluation result should be a rvalue
       call set_rvalue(res)
    endif

    !call disp_info(A, 'A3 = ')    
    
    if(associated(res) &
         .and. res%data == 0) then
       call array_duplicate(res, ops_array(1)%ptr)
    endif

    ! call push_back(buf, A%args(1))
    ! call push_back(buf_list, buf)
    ! buf_ptr => buf_data(buf_list)

    ! print*, "A%node_type = ", A%node_type
    kernel_addr = op_desc_list(A%node_type)%func
    ! write(*, "(A, Z16.16)"), "kernel_addr = ", kernel_addr
    
    num_oprands = size(ops_array)

    ! call disp_info(ops_array(1)%ptr, 'ops_array(1)%ptr = ')
    ! call disp_info(ops_array(2)%ptr, 'ops_array(2)%ptr = ')

    write(*, '(A, A, A, Z16.16)') "invoke kernel of ", &
         trim(op_desc_list(A%node_type)%kernel), &
         "  Addr : ", kernel_addr
    
    call invoke_kernel(kernel_addr, loc(res), &
         loc(ops_array), loc(ops_alpha), &
         loc(ops_beta), loc(A%args), &
         loc(num_oprands), loc(A%num_args), loc(ierr))

    ! call disp(res, 'res = ')
    print*, "evalation of ", trim(op_names(A%node_type)), " is complete."
    
    ! do i = 1, size(ops_array)
    !     call release_ptr(ops_array(i)%ptr)
    ! enddo
       
  end subroutine
  
  !> this function will be called in the following two situations:
  !> 1) A = B, where A and B are both arrays
  !> 2) func(A), where A is a array and passed to a function as an argument.
  !> in both case, we do not need a deep copy
  subroutine array_assign_array(A, B)
    implicit none
    type(array), intent(inout) :: A
    type(array), intent(in)  :: B
    integer :: ierr
    
    ! write(*,'(A, I0.1, A)'), '[', &
    !      get_rank(MPI_COMM_WORLD), '] &
    !      & called array_assign_array'

    if(is_rvalue(B)) then
       call array_copy(A, B)
    else
       call array_deep_copy(A, B)
    endif
  end subroutine

#:for type in [['i4', 'integer'],['r4', 'real'],['r8', 'real(8)']]
  subroutine array_assign_${type[0]}$(A, B)
    implicit none
    type(array), intent(inout) :: A
    type(${type[1]}$), intent(in) :: B
    integer :: ierr
    call set(A, B)
  end subroutine
#:endfor
  
  !> get number of sub nodes
  recursive function expr_sub_cnt(o, arg_mode, arg_is_root) result(res)
    implicit none
    type(node), intent(in) :: o
    integer, save :: cnt = 0
    integer :: i, res
    integer, optional :: arg_mode
    integer :: mode
    logical, optional :: arg_is_root
    logical :: is_root

    is_root = .true.
    if(present(arg_is_root)) is_root = arg_is_root
    
    res = 0    
    mode = 0
    if(present(arg_mode)) mode = arg_mode

    if(is_data(o)) then
       return
    endif

    if(mode == 1) then
       if(.not. is_ew(o)) then
          return
       endif
    endif

    do i = 1, size(o%operands)
       res = res + expr_sub_cnt(o%operands(i)%ptr, mode, .false.)
    end do
    res = res + size(o%operands)

    if(is_root) res = res + 1

  end function

  !> get the tree depth, skiping the non-ew node
  recursive function get_depth(o) result(res)
    implicit none
    type(node), intent(in) :: o
    integer :: res
    integer, allocatable :: sub_depth(:)
    integer :: i
    
    if(is_data(o) .or. (.not. is_ew(o))) then
       res = 1
       return
    endif

    allocate(sub_depth(size(o%operands)))
    do i = 1, size(o%operands)
       sub_depth(i) = get_depth(o%operands(i)%ptr)
    end do
    res = maxval(sub_depth) + 1
  end function

  !> general kernel infomation
  recursive subroutine gen_kernels(o, arg_is_root)
    implicit none
    type(node), intent(in) :: o
    character(len=MAX_TREESTR_LENGTH) :: str
    logical, optional :: arg_is_root
    logical :: is_root
    integer ::out_unit = 10
    integer :: i

    if(present(arg_is_root)) then
       is_root = arg_is_root
    else
       is_root = .true.
    endif

    if(get_rank(MPI_COMM_WORLD) == 0) then
       if(is_data(o)) return

       if(.not. is_ew(o)) then
          do i = 1, size(o%operands)
             call gen_kernels(o%operands(i)%ptr, .true.)
          enddo
          return
       endif

       if(is_root) then
          if(get_depth(o) >= 3) then
             call tree_to_string(str, o, SKIP_EW_NODE)
             ! print*, "node_type = ", op_names(o%node_type)
             ! print*, "depth = ", get_depth(o)
             ! print*, "str = ", trim(str)
             ! print*, "hash = ", djb_hash(trim(str))

             open(unit=out_unit,file=kernel_file, &
                  action='write', status='unknown', &
                  position='append')

             write(out_unit, "(I0.1,A)") &
                  abs(djb_hash(trim(str))), ":"//trim(str)

             close(out_unit)
          endif

          do i = 1, size(o%operands)
             call gen_kernels(o%operands(i)%ptr, .false.)
          enddo
       endif
    endif
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

       ! call reset_global_id()
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
#:for type1 in ['real(8)', 'real', 'integer', 'array', 'node']
#:for type2 in ['real(8)', 'real', 'integer', 'array', 'node']
#:for op in ['plus', 'minus', 'mult', 'divd', 'gt', 'ge', 'lt', 'le', 'eq', 'ne']
#:if any(type1 in s for s in ['real(8)', 'real', 'integer']) & 
     and any(type2 in s for s in ['real(8)', 'real', 'integer'])  
#:else  
#:set name1 = re.sub('\(8\)', '8', type1)
#:set name2 = re.sub('\(8\)', '8', type2)
  function expr_${name1}$_${op}$_${name2}$ (A, B) result(res)
    implicit none       
    type(${type1}$),intent(in)  :: A
    type(${type2}$),intent(in)  :: B
    type(node), pointer :: res
    type(node), pointer :: C, D

    
#:if type1 == 'node'
    call bind_ptr(C, A)
#:else
    allocate(C)
    call node_new(C, A)
#:endif

#:if type2 == 'node'
    call bind_ptr(D, B)
#:else
    allocate(D)
    call node_new(D, B)
#:endif

#:if type1 == 'node' and type2 == 'node'
    if(A%node_type == type_pow .and. B%node_type == type_scalar) then
       call disp_info(C, 'C = ')
       call disp_info(D, 'D = ')       
    endif
#:endif
    
    allocate(res)
    call node_new(res, type_${op}$, C, D)

#:if type1 == 'node' and type2 == 'node'    
    if(A%node_type == type_pow .and. B%node_type == type_scalar) then
       call disp_info(C, 'C = ')
       call disp_info(D, 'D = ')
       call disp_info(res, 'res = ')
    endif
#:endif
    
    call set_rvalue(res)
  end function
  
#:endif
#:endfor
#:endfor
#:endfor

#:for e in L
#:if e[1] >= 111
  #:for type in ['array', 'node']
  !> function ${e[2]}$(${type}$), return a tree node 
  function expr_${e[2]}$_${type}$ (o) result(res)
    implicit none
    type(${type}$), intent(in) :: o
    type(node), pointer :: res
    type(node), pointer :: tmp_o

#:if type == 'node'
    call bind_ptr(tmp_o, o)
#:else
    allocate(tmp_o)
    call node_new(tmp_o, o)
#:endif

    allocate(res)
    call node_new(res, ${e[0]}$, tmp_o)

    ! print*, "calling expr_${e[2]}$_${type}$ ..."
    ! print*, op_names(res%node_type),"==>", op_names(tmp_o%node_type)
    ! if(tmp_o%node_type == type_data) &
    !      print*, "xxxxxxxxxxxxxx = ", tmp_o%ref_cnt
    
    call set_rvalue(res)
  end function

#:endfor
#:endif
#:endfor

#:for type1 in ['node', 'array']
#:for type2 in ['integer', 'real', 'real(8)']
#:set name2 = re.sub('\(8\)', '8', type2)
  !> power function ${type1}$**${type2}$
  function expr_${type1}$_pow_${name2}$(o, n) result(res)
    implicit none
    type(${type1}$), intent(in) :: o
    type(${type2}$), intent(in) :: n
    type(node), pointer :: res
    type(node), pointer :: tmp_o

#:if type1 == 'node'
    call bind_ptr(tmp_o, o)
#:else
    allocate(tmp_o)
    call node_new(tmp_o, o)
#:endif

    allocate(res)
    call node_new(res, type_pow, tmp_o)
    res%num_args = 1
    res%args(1) = n
    call set_rvalue(res)
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
       
#:for data in ['node', 'array']  
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
    type(${data}$), intent(in) :: obj    
    ${ta[1]}$ :: a
    ${tb[1]}$ :: b    
    ${tc[1]}$ :: c
    type(node), pointer :: res, tmp
    integer :: dim
    integer :: xs, xe, ys, ye, zs, ze
    integer :: dim_x, dim_y, dim_z
    integer :: src_shape(3)

    allocate(res)

    ! call assert(is_valid(obj), &
    !      __FILE__, __LINE__, &
    !      'the input array/node must has a valid shape')
    
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
    
#:if data == 'array'
    !call node_new(res, obj) !create a data node
    call assign_ptr(res%data, obj)
    !call node_add_operand(res, tmp)
#:else
    call assign_ptr(tmp, obj) 
    call node_add_operand(res, tmp)
#:endif

    res%m_shape = (/dim_x, dim_y, dim_z/)
    res%node_type = type_ref
    call set_rvalue(res)
    
  end function
  
#:endfor
#:endfor
#:endfor
#:endfor  

#:for src_type in [['array', 'array', 'array'], &
       ['node', 'node', 'node'], &
       ['int', 'integer', 'scalar'], &
       ['real', 'real', 'scalar'], &
       ['real8', 'real(8)', 'scalar']]
#:for dst_type in [['array', 'array'], ['ref', 'node']]
#:set func_name ='set_{0}_from_{1}'.format(dst_type[0], src_type[0])
  subroutine ${func_name}$ (dst, src)
    type(${dst_type[1]}$), intent(in) :: dst
    type(${src_type[1]}$), intent(in) :: src
    type(${dst_type[1]}$), pointer :: dst_ptr
    type(${src_type[1]}$), pointer :: src_ptr
    type(ref_info) :: set_ref

    !print*, "calling function ${func_name}$"
    
#:if dst_type[0]=='array'  and src_type[2] == 'scalar'
    call range_to_ref (set_ref, range_all, range_all, range_all)
    call data_set_scalar (dst, real(src, 8), set_ref)
#:endif
    
#:if dst_type[0] == 'ref' and src_type[2] == 'scalar'
    call data_set_scalar (dst%data, real(src, 8),  dst%ref)
#:endif

#:if dst_type[0] == 'array' and src_type[2] == 'array'
    call bind_ptr(dst_ptr, dst)
    call bind_ptr(src_ptr, src)
    dst_ptr = src_ptr
#:endif

#:if dst_type[0] == 'array' and src_type[2] == 'node'
    call bind_ptr(dst_ptr, dst)
    call bind_ptr(src_ptr, src)
    dst_ptr = src_ptr
#:endif
    
#:if dst_type[0] == 'ref' and src_type[2] == 'array'
    call assert(is_ref(dst), __FILE__, __LINE__,&
         "the target must be a reference.")
    call set(dst, slice(src, range_all, range_all, range_all))
#:endif

#:if dst_type[0] == 'ref' and src_type[1] == 'node'
    call assert(associated(dst%data), __FILE__, __LINE__,&
         "as lvalue, node%data must be associated.")

    call assert(is_ref(dst), __FILE__, __LINE__,&
         "the target must be a reference.")

    call assign_ptr(src_ptr, src)
    
    if(allocated(src%operands) .and. &  !for expression node    
         size(src%operands)>1) then

       print*, "expr_node"

       !if the source node is an expression node,
       !evaluate it firt
       call eval(src_ptr%data, src_ptr, ierr)

       call data_set_ref_ref(dst%data, &
            src%data, dst%ref, src%ref)
       
    else if(src%node_type == type_data) then !for data node
       print*, "data_node"       
       !if the source node is an data node,
       !we have to create reference info using its shape
       call data_set_ref_ref(dst%data, &
            src_ptr%data, shape_to_ref(shape(src_ptr)), &
            src_ptr%ref)
    else
       ! call disp_info(dst, 'dst = ')
       ! call disp_info(src, 'src = ')
       ! call disp_info(dst%data, 'dst%data = ')
       ! call disp_info(src%data, 'src%data = ')

       call data_set_ref_ref(dst%data, &
            src_ptr%data, dst%ref, src_ptr%ref)
       
    end if
    !call disp_info(src_ptr, 'src_ptr = ')
    call release_ptr(src_ptr)
#:endif

    !print*, "finshing calling function ${func_name}$"
  end subroutine
#:endfor
#:endfor
  
  ! function replace_text (s,text,rep)  result(outs)
  !   character(*)        :: s,text,rep
  !   character(len(s)+100) :: outs     ! provide outs with extra 100 char len
  !   integer             :: i, nt, nr

  !   outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
  !   do
  !      i = index(outs,text(:nt)) ; if (i == 0) exit
  !      outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
  !   end do
  ! end function

  recursive subroutine gen_hash(o, arg_is_root)
    implicit none
    type(node), intent(inout) :: o
    integer :: mode
    character(len=MAX_TREESTR_LENGTH) :: str
    logical, optional :: arg_is_root
    logical :: is_root
    integer :: i
    
    if(present(arg_is_root)) then
       is_root = arg_is_root
    else
       is_root = .true.
    end if

    if(is_data(o)) return

    if((.not. is_ew(o))) then
       do i = 1, size(o%operands)
          call gen_hash(o%operands(i)%ptr, .true.)
       end do
       return
    endif
    
    if(is_root) then
       if(get_depth(o) >= 3) then
          call tree_to_string(str, o, SKIP_EW_NODE)
          o%hash = abs(djb_hash(trim(str)))
          ! print*, 'o%hash = ', o%hash
          ! print*, 'str = ', trim(str)
       endif
       
       do i = 1, size(o%operands)
          call gen_hash(o%operands(i)%ptr, .false.)
       enddo
    endif
  end subroutine
  
  recursive subroutine tree_to_string(str, o, arg_mode)
    implicit none
    character(len=MAX_TREESTR_LENGTH), intent(out) :: str
    type(node), intent(in) :: o
    character(len=10) :: expr_form
    integer :: i, num_operands, i1,i2
    character(len=MAX_TREESTR_LENGTH), allocatable :: ops_strs(:)
    integer, optional :: arg_mode
    integer :: mode

    if(present(arg_mode)) then
       mode = arg_mode
    else
       mode = ALL_NODE
    endif

    if(is_data(o)) then
       str = '(X*A+Y)'
       return
    endif

    if(mode == SKIP_EW_NODE .and. &
         (.not. is_ew(o))) then
       str = '(X*A+Y)'
       return
    endif

    num_operands = size(o%operands)
    
    allocate(ops_strs(num_operands))
    
    do i = 1, size(o%operands)
       call tree_to_string(ops_strs(i), &
            o%operands(i)%ptr, mode)
    enddo
    
    expr_form = op_desc_list(o%node_type)%expr_form

    select case (num_operands)
    case (1)
       i1 = index(expr_form, 'A')
       ! print*, "expr_form = ", expr_form, "i = ", i1
       ! print*, "expr_form(:i1-1) = ", expr_form(:i1-1)
       ! print*, "expr_form(i1+1:) = ", expr_form(i1+1:)
       write(str, "(A, A, A)") &
            "X*(", expr_form(:i1-1)//trim(ops_strs(1))&
            //trim(expr_form(i1+1:)), ")+Y"

    case (2)
       i1 = index(expr_form, 'A')
       i2 = index(expr_form, 'B')
       write(str, "(A,A,A)") "X*(", &
            expr_form(:i1-1)//trim(ops_strs(1))&
            //expr_form(i1+1:i2-1) &
            //trim(adjustl(ops_strs(2)))//trim(expr_form(i2+1:)), ")+Y"
    end select
  end subroutine

  
end module

