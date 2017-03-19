module ot_expr
  use ot_common
  use ot_tensor
  use dispmodule
  use ot_ref
  use ot_vector
  use ot_geom
  use ot_node
  
#include "node_type.h"
!#define ASSERT(x,msg) call assert(x, __FILE__, __LINE__, msg)
#:set DEBUG = 1
#define DEBUG
#:include "type_def.fypp"
  
  public :: operator(+),operator(-)
  public :: operator(*),operator(/)
  public :: assignment(=)

  interface display
     module procedure display1, display2
  end interface
  
  interface assignment (=)
     module procedure tensor_assign_tensor
     module procedure node_assign_tensor
     module procedure node_assign_node
  end interface assignment (=)


  interface slice
#:for data in ['node', 'tensor']
#:for ta in [['int', 'integer'], ['arr', 'integer, dimension(:)']]
#:for tb in [['int', 'integer'], ['arr', 'integer, dimension(:)']]
#:for tc in [['int', 'integer'], ['arr', 'integer, dimension(:)']]    
     module procedure slice_${data}$_${ta[0]}$_${tb[0]}$_${tc[0]}$
#:endfor
#:endfor
#:endfor
#:endfor
  end interface slice

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

  
contains

  subroutine init_expr(ierr)
    use ot_tensor
    implicit none
    integer, intent(out) :: ierr
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
    A%m_dim = B%m_dim
    A%m_shape = B%m_shape
    A%alpha = B%alpha
    A%beta = B%beta
    A%scalar = B%scalar
    A%is_implicit = B%is_implicit
    call copy_node_ptr(A%operands, B%operands)
    call copy_node_ptr(A%opt_operands, B%opt_operands)
    
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

    call node_combine(C, res)
    
    call eval(A, C, ierr)

    !call display(C%data)
    
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

  subroutine set(obj1, obj2)
    implicit none
    type(node), intent(inout) :: obj1
    type(node), intent(in) :: obj2
    ! ref_node, ref_node
    ! node, ref_node
    ! ref_node, node
    ! node, node
  end subroutine

  !*****************************
  !> slice functions
  !*****************************
#:set itype=[['int', 'integer'], &
       ['range', 'type(range)'], &
       ['arr', 'integer,dimension(:)']]
       
#:for data in ['node', 'tensor']  
#:for ta in itype
#:for tb in itype
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

    allocate(res)
    
#:if ta[0] == 'int'
    xs = a; xe = a;
    dim_x = xe - xs + 1
#:elif ta[0] == 'range'
    xs = a%lower; xe = a%upper
    dim_x = xe - xs + 1
#:else
    allocate(res%ref%ix(size(a)))
    res%ref%ix = a
    dim_x = size(a)
#:endif

#:if tb[0] == 'int'
    ys = b; ye = b;
    dim_y = ye - ys + 1
#:elif tb[0] == 'range'
    ys = b%lower; ye = b%upper
    dim_y = ye - ys + 1    
#:else
    allocate(res%ref%iy(size(b)))
    res%ref%iy = b
    dim_y = size(b)
#:endif

#:if tc[0] == 'int'
    zs = c; ze = c;
    dim_z = ze - zs + 1    
#:elif tc[0] == 'range'
    zs = c%lower; ze = c%upper
    dim_z = ze - zs + 1    
#:else
    allocate(res%ref%iz(size(c)))
    res%ref%iz = c
    dim_z = size(c)
#:endif
    
#:if data == 'tensor'
    allocate(tmp)
    call node_new(tmp, obj) !create a data node
    call node_add_operand(res, tmp)
#:else
    call assign_ptr(tmp, obj) 
    call node_add_operand(res, tmp)
#:endif

    res%m_dim = 3
    res%m_shape = (/dim_x, dim_y, dim_z/)
    res%node_type = type_ref

  end function
  
#:endfor
#:endfor
#:endfor
#:endfor  

  !*********************************
  !> evaluate expression nodes
  !*********************************
  recursive subroutine eval(res, A, ierr, is_root_arg)
    implicit none
    type(tensor), intent(inout) :: res
    type(node), intent(in) :: A
    integer :: i, num_oprands
    integer, intent(out) :: ierr
    real(8), allocatable :: ops_alpha(:), ops_beta(:)
    real(8) :: alpha, beta
    type(tensor_ptr), allocatable :: ops_tensor(:)
    logical, optional, intent(in) :: is_root_arg
    logical :: is_root
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

    if(present(is_root_arg)) then
       is_root = is_root_arg
    else
       is_root = .true.
    endif

    call disp_info(A, 'A = ')
    
    !if this node is the root node, force to evaluate
    if(is_data(A) .and. .not. is_root) return

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

    if(is_root) then
       alpha = A%alpha
       beta = A%beta
    else
       alpha = 0
       beta = 1.0
    endif
    
    if(need_eval(A) .or. is_root) then
       select case(A%node_type)
#:for e in L
#:if e[4]=='A' or e[4]=='B' or e[4]=='C'
       case (${e[1]}$)
          call ${e[2]}$_tensors(res, &
               ops_tensor, ops_alpha, ops_beta, A%args, alpha, beta)
#:endif
#:endfor
       case default
          call eval_tensors(res, &
               ops_tensor, ops_alpha, ops_beta, A%args, alpha, beta)
       end select

       ! if(.not. associated(A%data)) then
       !    print*, "execute function ", op_names(A%node_type), " failed."
       !    call abort()
       ! end if
    endif
      
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
    
    if(B%is_implicit) then
       call tensor_copy(A, B)
    else
       call tensor_deep_copy(A, B)
    end if
  end subroutine
  
  recursive subroutine node_combine(A, res) 
    implicit none
    type(node), intent(inout), target :: A
    type(node_ptr),  allocatable, intent(out) :: res(:)
    type(node_ptr),  allocatable :: left(:), right(:)
    type(node_ptr),  allocatable :: subs(:), tmp(:)
    integer i, cnt

    cnt = 0
    
    if(A%node_type == type_data &
         .or. A%node_type==type_scalar) return

    do i = 1,size(A%operands)
       !combine the left node
       call node_combine(A%operands(i)%ptr, subs)

       ! if the node_type does not match its leaf node_type,
       ! not accept the returned list from leaf
       if(A%node_type /= A%operands(i)%ptr%node_type .or. &
            (.not. is_arithmetic(A))) then
          if(allocated(subs)) deallocate(subs)
          allocate(subs(1))
          call assign_ptr(subs(1)%ptr, A%operands(i)%ptr)
       end if
       
       if(allocated(subs)) then
          cnt = cnt + size(subs)

          call reallocate(res, cnt)
          res(cnt-size(subs)+1 : cnt) = subs
       endif
    enddo

    call copy_node_ptr(A%opt_operands, res)

! #ifdef DEBUG
!     write(*,*) "**********************"
!     write(*, "(1X, Z16.16, A, I0.3)"), loc(A), " type=", A%node_type
    
!     if(allocated(A%opt_operands)) then
!        do i = 1, size(A%opt_operands)
!           write(*, "(5X, Z16.16, A, I0.3)"), &
!                loc(A%opt_operands(i)%ptr), &
!                " type=", A%opt_operands(i)%ptr%node_type 

!        enddo
!     endif
! #endif
    
  end subroutine


!   recursive subroutine node_optimize(A, res) 
!     implicit none
!     type(node), intent(inout), target :: A
!     type(node_ptr),  allocatable, intent(out) :: res(:)
!     type(node_ptr),  allocatable :: subs(:), tmp(:)
!     integer i, cnt, num_oprands
!     type(node), pointer :: child
    
!     cnt = 0
    
!     if(is_data(A)) return

!     num_oprands = size(A%operands)

!     do i = 1, num_oprands

!        child => A%operands(i)%ptr
       
!        if(A%node_type == child%node_type .and. &
!             is_arithmetic(A)) then
!           call reallocate(A%operands, num_oprands+size(A%operands))
!           A%operands()
!        end if
       
!        ! if the node_type does not match its leaf node_type,
!        ! not accept the returned list from leaf
!        if(A%node_type == A%operands(i)%ptr%node_type) then
!           if(allocated(subs)) deallocate(subs)
!           allocate(subs(1))
!           subs(1)%ptr => A%operands(i)%ptr
!        end if
       
!        if(allocated(subs)) then
!           cnt = cnt + size(subs)

!           call reallocate(res, cnt)
!           res(cnt-size(subs)+1 : cnt) = subs
!        endif
!     enddo

!     do i = 1, size(A%operands) 
!        !combine the left node
!        call node_optimize(A%operands(i)%ptr, subs)
!     enddo
!     call copy_node_ptr(A%opt_operands, res)

! #ifdef DEBUG
!     write(*,*) "**********************"
!     write(*, "(1X, Z16.16, A, I0.3)"), loc(A), " type=", A%node_type
    
!     if(allocated(A%opt_operands)) then
!        do i = 1, size(A%opt_operands)
!           write(*, "(5X, Z16.16, A, I0.3)"), &
!                loc(A%opt_operands(i)%ptr), &
!                " type=", A%opt_operands(i)%ptr%node_type 

!        enddo
!     endif
! #endif
!    end subroutine
  

#:for o in ['', 'opt_']
  recursive subroutine write_${o}$graph(A, is_root, file)
    implicit none
    type(node), intent(inout) :: A
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

    id = get_global_id()
    A%id = id
    
#:set SHOW_ALPHA_BETA = 1
    
#:if SHOW_ALPHA_BETA > 0
    write(out_unit, format_label) id, &
         '[label="[',trim(adjustl(op_name)),']\n id=', A%id, '\n(',&
         A%alpha, ',',A%beta, ')"];'
#:else
    write(out_unit, format_label1) id, &
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
    
! #ifdef DEBUG
! #:if o == "opt_"
!     print*, "-------------------------"
! #:else
!     print*, "~~~~~~~~~~~~~~~~~~~~~~~~~"
! #:endif
!     do i = 1, size(A%${o}$operands)
!        write(*, "(I0.1, A, Z16.16)"), &
!             i," = ", loc(A%${o}$operands(i)%ptr)           
!     end do
! #endif

    if(allocated(A%${o}$operands)) then
       do i = 1, size(A%${o}$operands)
          call write_${o}$graph(A%${o}$operands(i)%ptr,  .false.)
          write(out_unit, format_map), &
               A%id, "->", A%${o}$operands(i)%ptr%id, ";"
       end do
    endif
    
    if(am_i_root) then
       write(out_unit, *) "}"
       close (out_unit)
    end if
    
  end subroutine
#:endfor

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
end module

