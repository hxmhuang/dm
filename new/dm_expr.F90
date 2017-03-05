module dm_expr
  use dm_common
  use dm_tensor
#include "dm_type.h"
!#define ASSERT(x,msg) call assert(x, __FILE__, __LINE__, msg)
#:set DEBUG = 1
  
  public :: operator(+),operator(-)
  public :: operator(*),operator(/)
  public :: assignment(=)

  type node_array
     type(node), pointer :: ptr => null()
  end type node_array
  
  type node
     type(tensor), pointer :: data => null()

     !the operands for a operation, e.g. A + B,, A and B are operands
     !for a binary operator, the number operands is always 2
     type(node_array), allocatable, dimension(:) :: operands

     !optimized operands, e.g. A + B + C, opt_operands are (A,B,C)
     type(node_array), allocatable, dimension(:) :: opt_operands
     
     !node type, 0 for data, 1 for '+', 2 for '-',
     !3 for '*', 4 for '/', etc.
     integer :: node_type

     ! the node shape inherit from its left/right
     ! expression or the tensor within it
     integer :: m_dim = 0
     integer :: m_shape(3) = (/0,0,0/)
     
     real(8) :: alpha = 0
     real(8) :: beta  = 1
     real(8) :: scalar = 0
     logical :: is_implicit
     integer :: id
  end type node

  !the following code using preprossor to create interfaces
#:for op in [['plus','+'], ['minus','-'], ['mult','*'], ['divd','/']]
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
  
  interface display
     module procedure display1, display2
  end interface
  
  interface assignment (=)
     module procedure tensor_assign_tensor
     module procedure node_assign_tensor
     module procedure node_assign_node
  end interface assignment (=)

  interface node_new
     module procedure node_new1, node_new2, node_new3
     module procedure node_new4, node_new5, node_new6
  end interface node_new
  
  interface reallocate
     module procedure reallocate1,reallocate2
  end interface reallocate

  interface operator (-)
     module procedure uminus_node, uminus_tensor
  end interface operator (-)

  interface operator (+)
     module procedure uplus_node, uplus_tensor
  end interface operator (+)
  
contains

  !the following code using preprossor to create subroutines  
#:for type1 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:for type2 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:for op in ['plus', 'minus', 'mult', 'divd']
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
#:if (name2 == 'real8')
    write(*,*) "AA=", C%node_type
    write(*,*) "C=",  C%m_shape
#:endif
    call node_new(res, type_${op}$, C, D)
  end function
  
#:endif
#:endfor
#:endfor
#:endfor

  function uminus_node (A) result(res)
    implicit none       
    type(node),intent(in),target  :: A
    type(node), target :: res
    type(node), pointer :: C, D

    allocate(C, D)
    
    call node_new(C, 0)
    call node_new(D, A)
    
    call node_new(res, type_minus, C, D)
  end function

  function uminus_tensor (A) result(res)
    implicit none       
    type(tensor),intent(in),target  :: A
    type(node), target :: res
    type(node), pointer :: C, D

    allocate(C, D)
    
    call node_new(C, 0)
    call node_new(D, A)
    
    call node_new(res, type_minus, C, D)
  end function

  function uplus_node (A) result(res)
    implicit none       
    type(node),intent(in),target  :: A
    type(node), target :: res
    type(node), pointer :: C, D

    allocate(C, D)
    
    call node_new(C, 0)
    call node_new(D, A)
    
    call node_new(res, type_plus, C, D)
  end function

  function uplus_tensor (A) result(res)
    implicit none       
    type(tensor),intent(in),target  :: A
    type(node), target :: res
    type(node), pointer :: C, D

    allocate(C, D)
    
    call node_new(C, 0)
    call node_new(D, A)
    
    call node_new(res, type_plus, C, D)
  end function
  
  subroutine node_new1(B, A) 
    implicit none
    type(tensor), intent(in), target :: A
    type(node),   intent(out) :: B
    B%m_dim = A%m_dim
    B%m_shape = A%m_shape
    B%data => A
    B%node_type = type_data
  end subroutine
  
  subroutine node_new2(B, A) 
    implicit none
    type(node), intent(in)   :: A
    type(node), intent(out)  :: B
    integer i
    
    B%m_dim = A%m_dim
    B%m_shape = A%m_shape
    B%alpha = A%alpha
    B%beta = A%beta
    B%data => A%data

    B%scalar = B%scalar
    B%node_type = A%node_type

    call copy_node_array(B%opt_operands, A%opt_operands)
    call copy_node_array(B%operands, A%operands)
  end subroutine

  subroutine node_new3(C, op_type, left, right)
    implicit none
    integer, intent(in) :: op_type
    type(node), intent(in), target :: left, right
    type(node), intent(out) :: C
    real(8) :: x

    if((.not. is_scalar(left)) .and. (.not. is_scalar(right))) then
       call assert(left%m_dim == right%m_dim .and. &
            left%m_dim == right%m_dim, &
            __FILE__, __LINE__, &
            "shape of the left and right leaf does not match.")
    endif
    
    if(is_scalar(left)) then
       !copy the right to C
       C%m_dim   = right%m_dim
       C%m_shape = right%m_shape
       
       x = left % scalar
       select case(op_type)
       case (type_plus) ! x + A
          C%alpha = x + right%alpha
          C%beta  = right%beta
          C%node_type = right%node_type
          C%data => right%data
          C%is_implicit = right%is_implicit          
          call copy_node_array(C%operands, right%operands)
       case (type_minus) ! x - A
          C%alpha = x - right%alpha
          C%beta  = -right%beta
          C%node_type = right%node_type
          C%data => right%data
          C%is_implicit = right%is_implicit                    
          call copy_node_array(C%operands, right%operands)          
       case (type_mult) ! x .* A
          C%alpha = x * right%alpha
          C%beta  = x * right%beta
          C%node_type = right%node_type
          C%data  => right%data
          C%is_implicit = right%is_implicit                    
          call copy_node_array(C%operands, right%operands)
       case (type_divd) ! x ./ A
          C%beta = x
          C%node_type = type_rcp
          allocate(C%operands(1))
          C%operands(1)%ptr => right
          C%is_implicit = .true.
       end select
       
    else if(is_scalar(right)) then
       x = right % scalar
       select case(op_type)
       case (type_plus) ! A + x
          C%alpha = left%alpha + x
          C%beta  = left%beta
          C%node_type = left%node_type
          call copy_node_array(C%operands, left%operands)
       case (type_minus) ! A - a
          C%alpha = left%alpha - x
          C%beta  = left%beta
          C%node_type = left%node_type
          call copy_node_array(C%operands, left%operands)
       case (type_mult) ! A .* a
          C%alpha = left%alpha * x
          C%beta  = left%beta  * x
          C%node_type = left%node_type
          call copy_node_array(C%operands, left%operands)
       case (type_divd)
          C%alpha = left%alpha / x
          C%alpha = left%beta  / x
          C%node_type = left%node_type
          call copy_node_array(C%operands, left%operands)
       end select

       !copy the left to C
       C%m_dim   = left%m_dim
       C%m_shape = left%m_shape
       C%is_implicit = left%is_implicit
       call copy_node_array(C%opt_operands, left%opt_operands)
    else
       allocate(C%operands(2))
       C%operands(1)%ptr => left
       C%operands(2)%ptr => right
       
       C%m_dim   = left%m_dim
       C%m_shape = left%m_shape
       C%node_type = op_type
       C%is_implicit = .true.
    end if
  end subroutine

  subroutine node_new4(B, scalar) 
    implicit none
    real(8), intent(in), target :: scalar
    type(node), intent(out) :: B
    
    B%m_dim = 0
    B%m_shape = (/1, 1, 1/)
    B%scalar = scalar
    B%node_type = type_scalar
  end subroutine

  subroutine node_new5(B, scalar)
    implicit none
    real, intent(in):: scalar
    type(node),intent(out) :: B

    call node_new4(B, real(scalar,8))
  end subroutine

  subroutine node_new6(B, scalar)
    implicit none
    integer, intent(in) :: scalar
    type(node), intent(out) :: B

    call node_new4(B, real(scalar,8))
  end subroutine
  
  subroutine display1(A, msg)
    implicit none
    type(node) :: A
    character(len=*),intent(in),optional :: msg
    if(present(msg)) then
       call display2(A%data, msg)
    else
       call display2(A%data)
    end if
  end subroutine

  ! recursive subroutine eval(B)
  !   implicit none
  !   type(node), intent(inout) :: B

  !   if(B%node_type == type_data) return 

  !   call assert(associated(B%left) .and. &
  !        associated(B%right), __FILE__, __LINE__, &
  !        "pointer not associated.")

  !   !evaluate left and right leaves recursively.
  !   if(B%left%node_type  /= type_data) call eval(B%left)
  !   if(B%right%node_type /= type_data) call eval(B%right)       

  !   select case(B%node_type)
  !   case (type_plus)
  !      !C = A + B       
  !      call tensor_plus(B%data, B%left%data, B%right%data)  
  !   case (type_minus)
  !      !C = A - B       
  !      call tensor_minus(B%data, B%left%data, B%right%data) 
  !   case (type_mult)
  !      !C = A .* B       
  !      call tensor_mult(B%data, B%left%data, B%right%data)  
  !   case (type_divd)
  !      !C = A ./ B       
  !      call tensor_divd(B%data, B%left%data, B%right%data)  
  !   end select
    
  !   B%node_type = type_data
  ! end subroutine 

  subroutine node_assign_node(A, B)
    implicit none
    type(node), intent(inout), target :: A
    type(node), intent(in),    target :: B

#:if DEBUG > 0
    print*, "call node_assign_node"
#:endif
    
    A%data  => B%data
    A%node_type = B%node_type
    A%m_dim = B%m_dim
    A%m_shape = B%m_shape
    A%alpha = B%alpha
    A%beta = B%beta
    A%scalar = B%scalar
    
    call copy_node_array(A%operands, B%operands)
    call copy_node_array(A%opt_operands, B%opt_operands)
    
  end subroutine
  
  subroutine node_assign_tensor(A, B)
    implicit none
    type(tensor), intent(out) :: A
    
    !the type of B must be "intent(in)"
    type(node), intent(in)   :: B 
    type(node) :: C
    type(node_array), allocatable :: res(:)
    integer ierr
    
    call node_new(C, B)

    call node_combine(C, res)
    
    call write_graph(C, .true., "graph.dot")
    call write_opt_graph(C, .true., "opt_graph.dot")
    
    ! call node_destroy(C, ierr)
    
    !write(*, "(I3)") size(C%opt_operands)
    !call eval(B)
    
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

  !> this function will be called in the following two situations:
  !> 1) A = B, where A and B are both tensors
  !> 2) func(A), where A is a tensor and passed to a function as an argument.
  !> in both case, we do not need a deep copy
  subroutine tensor_assign_tensor(A, B)
    implicit none
    type(tensor), intent(out) :: A
    type(tensor), intent(in) :: B

    if(B%is_implicit) then
       call tensor_copy(A, B)
    else
       call tensor_deep_copy(A, B)
    end if
  end subroutine
  
  recursive subroutine node_combine(A, res) 
    implicit none
    type(node), intent(inout), target :: A
    type(node_array),  allocatable, intent(out) :: res(:)
    type(node_array),  allocatable :: left(:), right(:)
    type(node_array),  allocatable :: subs(:), tmp(:)
    integer i, cnt

    cnt = 0
    
    if(A%node_type == type_data &
         .or. A%node_type==type_scalar) return

    do i = 1,size(A%operands)
       !combine the left node
       call node_combine(A%operands(i)%ptr, subs)

       !if the node_type does not match its leaf node_type,
       ! not accept the returned list from leaf
       if(A%node_type /= A%operands(i)%ptr%node_type) then
          if(allocated(subs))  deallocate(subs)
          allocate(subs(1))
          subs(1)%ptr => A%operands(i)%ptr
       end if
       
       if(allocated(subs)) then
          cnt = cnt + size(subs)

          call reallocate(res, cnt)
          res(cnt-size(subs)+1 : cnt) = subs
       endif
    enddo

    call copy_node_array(A%opt_operands, res)

#ifdef DEBUG
    write(*,*) "**********************"
    write(*, "(1X, Z16.16, A, I0.3)"), loc(A), " type=", A%node_type
    
    if(allocated(A%opt_operands)) then
       do i = 1, size(A%opt_operands)
          write(*, "(5X, Z16.16, A, I0.3)"), &
               loc(A%opt_operands(i)%ptr), &
               " type=", A%opt_operands(i)%ptr%node_type 

       enddo
    endif
#endif
    
  end subroutine

  recursive subroutine node_destroy(A, ierr)
    implicit none
    type(node), intent(inout), target :: A
    integer :: ierr
    integer i

    if(A%node_type == type_data .or. &
         A%node_type == type_scalar) return
    
    ! do i = 1, size(A%opt_operands)
    !    write(*, "(Z16.16, A, I4, A, I4)"), &
    !         loc(A%opt_operands(i)%ptr), " = ", A%opt_operands(i)%ptr%node_type, &
    !         " = ", size(A%opt_operands(i)%ptr%opt_operands)
    ! enddo
    
    if(allocated(A%operands)) then
       do i = 1, size(A%operands)
          call node_destroy(A%operands(i)%ptr, ierr)
       enddo
    endif
    
    !destroy the operands
    deallocate(A%opt_operands)
    !A%opt_operands => null()

    !destroy data if the data is implicit
    if(associated(A%data)) then
       if(A%data%is_implicit) then
          call tensor_destroy(A%data, ierr)
          A%data => null()
       end if
    end if
  end subroutine

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
         format_label = "(I0.1,A,A,A,F4.1,A,F4.1,A)"
    character(len=*), parameter :: &
         format_label1 = "(I0.1,A,A,A)"
    character(len=*), parameter :: &
         format_map = "(I0.1,A,I0.1,A)"
    
    id = get_global_counter()
    A%id = id
    
    if(present(is_root)) then
       am_i_root = is_root
    else
       am_i_root = .true.
    endif

    if(am_i_root) then
       if(present(file)) then
          open(unit=out_unit,file=file, &
               action="write", status="replace")
       else
          open(unit=out_unit,file="graph.dot", &
               action="write", status="replace")
       endif
       write(out_unit, *) "digraph G {"
    end if

    select case (A%node_type)
    case (type_scalar)
       write(op_name, "(F4.1)") A%scalar
    case (type_data)
       write(op_name, *) "Matrix"
    case (type_plus)
       write(op_name, *) "(+)"
    case (type_minus)
       write(op_name, *) "(-)"
    case (type_mult)
       write(op_name, *) "(.*)"
    case (type_divd)
       write(op_name, *) "(./)"
    case default
       write(op_name, *) "func"
    end select

#:set SHOW_ALPHA_BETA = 1
    
#:if SHOW_ALPHA_BETA > 0
    write(out_unit, format_label) id, &
         '[label="',trim(op_name), '\n(',&
         A%alpha, ',',A%beta, ')"];'
#:else
    write(out_unit, format_label1) id, &
         '[label="',trim(op_name), '"];'
#:endif
    
    !print*, label1
    if(A%node_type == type_data .or. A%node_type == type_scalar) then
       return
    end if
    
#ifdef DEBUG
#:if o == "opt_"
    print*, "-------------------------"
#:else
    print*, "~~~~~~~~~~~~~~~~~~~~~~~~~"
#:endif
    
    do i = 1, size(A%${o}$operands)
       write(*, "(I0.1, A, Z16.16)"), &
            i," = ", loc(A%${o}$operands(i)%ptr)           
    end do
#endif    
    do i = 1, size(A%${o}$operands)
       call write_${o}$graph(A%${o}$operands(i)%ptr,  .false.)
       write(out_unit, format_map), &
            A%id, "->", A%${o}$operands(i)%ptr%id, ";"
    end do
    
    if(am_i_root) then
       write(out_unit, *) "}"
       close (out_unit)
    end if

    
  end subroutine
#:endfor
  
  function is_scalar(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    res = (A%node_type == type_scalar)
  end function

  subroutine copy_node_array(dst, src)
    implicit none
    type(node_array), allocatable, intent(in) :: src(:)
    type(node_array), allocatable, intent(inout) :: dst(:)

    if(allocated(src)) then
       if(allocated(dst)) deallocate(dst)
       allocate(dst(size(src)))
       dst = src
    end if
  end subroutine

  !> reallocate an array while keeping its data
  subroutine reallocate1(arr, len)
    implicit none
    type(node_array), allocatable, intent(inout) :: arr(:)
    type(node_array), allocatable :: tmp(:)
    integer :: len, copy_len
    
    if(allocated(arr)) then
       !save the data to tmp
       allocate(tmp(size(arr)))
       tmp = arr
       deallocate(arr)
       
       !allocate new size
       allocate(arr(len))

       !copy data back from tmp
       copy_len = min(size(tmp), len)
       arr(1:copy_len) = tmp(1:copy_len)
       deallocate(tmp)
    else
       allocate(arr(len))
    endif
  end subroutine

  !> reallocate an array while keeping its data
  subroutine reallocate2(arr, len)
    implicit none
    type(node_array), pointer, intent(inout) :: arr(:)
    type(node_array), pointer :: tmp(:)
    integer :: len, copy_len
    
    if(associated(arr)) then
       !save the data to tmp
       allocate(tmp(size(arr)))
       tmp = arr
       deallocate(arr)
       
       !allocate new size
       allocate(arr(len))

       !copy data back from tmp
       copy_len = min(size(tmp), len)
       arr(1:copy_len) = tmp(1:copy_len)
       deallocate(tmp)
    else
       allocate(arr(len))
    endif
  end subroutine
  
end module

