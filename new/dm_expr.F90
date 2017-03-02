module dm_expr
  use dm_common
  use dm_tensor
#include "dm_type.h"
!#define ASSERT(x,msg) call assert(x, __FILE__, __LINE__, msg)
#define __NL__

  ! write(*, "(A, Z16.16)") "$1_$2_$3 : loc(A)=", loc(A)
  ! write(*, "(A, Z16.16)") "$1_$2_$3 : loc(B)=", loc(B)
  ! write(*, "(A, Z16.16)") "$1_$2_$3 : loc(C)=", loc(C(1))
  ! write(*, "(A, Z16.16)") "$1_$2_$3 : loc(D)=", loc(C(2))
  ! write(*, "(A, Z16.16)") "$1_$2_$3 : loc(res)=", loc(res)
  ! write(*, "(A, Z16.16)") "$1_$2_$3 : loc(res%left)=",  loc(res%left)
  ! write(*, "(A, Z16.16)") "$1_$2_$3 : loc(res%right)=", loc(res%right)
  
! !we use macro preprocessor to generate codes
! define(`evaluate', `eval')
! undefine(`eval')
! define(`CONCAT', `$1$2$3')dnl
! define(`OV', `
! define(`aa', ifelse($1, `real(8)', `real8', $1))dnl
! define(`bb', ifelse($3, `real(8)', `real8', $3))dnl
! function CONCAT(aa,_$2_,bb) (A, B) result(res)
!    implicit none       
!    type($1),intent(in),target  :: A
!    type($3),intent(in),target  :: B
!    type(node), target :: res
!    type(node), pointer :: C(:)

!    allocate(C(2))
   
!    C(1) = node_new(A)
!    C(2) = node_new(B)
!    res = node_new(type_$2, C(1), C(2))
! end function
! !')

! define(`OV_NAME',`module procedure CONCAT(ifelse($1, `real(8)', `real8', $1),_$2_,ifelse($3, `real(8)', `real8', $3))')

public :: operator(+),operator(-)
public :: operator(*),operator(/)
public :: assignment(=)

  type node_array
     type(node), pointer :: ptr => null()
  end type node_array
  
  type node
     type(tensor), pointer :: data => null()
     
     !the left and right leaf for the tree
     ! type(node), pointer :: left => null()
     ! type(node), pointer :: right => null()

     type(node_array), allocatable, dimension(:) :: subs

     type(node_array), allocatable, dimension(:) :: operands
     
     !node type, 0 for data, 1 for '+', 2 for '-',
     !3 for '*', 4 for '/'
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
    type(node), pointer :: C(:)

    allocate(C(2))

    call node_new(C(1), A)
    call node_new(C(2), B)
#:if (name2 == 'real8')
    write(*,*) "AA=", C(2)%node_type
    write(*,*) "C(2)=", C(2)%m_shape
#:endif
    call node_new(res, type_${op}$, C(1), C(2))
  end function
  
#:endif
#:endfor
#:endfor
#:endfor
  
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
    type(node), intent(in) :: A
    type(node), intent(out)  :: B
    integer i
    
    B%m_dim = A%m_dim
    B%m_shape = A%m_shape
    B%alpha = A%alpha
    B%beta = A%beta
    
    ! B%left => A%left
    ! B%right => A%right
    B%data => A%data

    if(B%node_type == type_scalar) print*, "node_new2 : B%scalar=", B%scalar
    B%scalar = B%scalar
    B%node_type = A%node_type

    !B%operands => A%operands
    ! if(associated(A%operands)) then
    !    allocate(B%operands(size(A%operands)))
    !    B%operands = A%operands
    ! end if

    call copy_node_array(B%operands, A%operands)
    call copy_node_array(B%subs, A%subs)
    
    ! if(associated(A%subs)) then
    !    allocate(B%subs(size(A%subs)))
    !    B%subs = A%subs

    !    ! do i = 1, size(A%subs)
    !    !    write(*, "(I0.1, A, Z16.16, A, I0.1)"), i," = ", loc(A%subs(i)%ptr), "type=",A%node_type           
    !    !call write_graph(A%subs(i)%ptr,  .false.)
    !    !write(out_unit, format_map), A%id, "->", loc(A%subs), ";"
    !    end do
       
    !end if
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
       x = left % scalar
       select case(op_type)
       case (type_plus) ! x + A
          C%alpha = x + right%alpha
          C%beta  = right%beta
          C%node_type = right%node_type
          C%data => right%data
          call copy_node_array(C%subs, right%subs)
       case (type_minus) ! x - A
          C%alpha = x - right%alpha
          C%beta  = -right%beta
          C%node_type = right%node_type
          C%data => right%data
          call copy_node_array(C%subs, right%subs)          
       case (type_mult) ! x .* A
          C%alpha = x * right%alpha
          C%beta  = x * right%beta
          C%node_type = right%node_type
          C%data  => right%data
          call copy_node_array(C%subs, right%subs)
       case (type_divd) ! x ./ A
          C%node_type = type_rcp
          allocate(C%subs(1))
          C%subs(1)%ptr => right
       end select

       !copy the right to C
       C%m_dim   = right%m_dim
       C%m_shape = right%m_shape
       C%is_implicit = right%is_implicit
       ! if(associated(right%operands)) then
       !    allocate(C%operands(size(right%operands)))
       !    C%operands = right%operands
       ! end if
    else if(is_scalar(right)) then
       x = right % scalar
       select case(op_type)
       case (type_plus) ! A + x
          C%alpha = left%alpha + x
          C%beta  = left%beta
          C%node_type = left%node_type
          
          call copy_node_array(C%subs, left%subs)
          
       case (type_minus) ! A - a
          C%alpha = left%alpha - x
          C%beta  = left%beta
          C%node_type = left%node_type
          
          call copy_node_array(C%subs, left%subs)

       case (type_mult) ! A .* a
          C%alpha = left%alpha * x
          C%beta  = left%beta  * x
          C%node_type = left%node_type
          
          call copy_node_array(C%subs, left%subs)

       case (type_divd)
          C%alpha = left%alpha / x
          C%alpha = left%beta  / x
          C%node_type = left%node_type
          
          call copy_node_array(C%subs, left%subs)
       end select

       !copy the left to C
       C%m_dim   = left%m_dim
       C%m_shape = left%m_shape
       C%data    => left%data
       C%is_implicit = left%is_implicit
       ! if(associated(left%operands)) then
       !    allocate(C%operands(size(left%operands)))
       !    C%operands = left%operands
       ! end if
       call copy_node_array(C%operands, left%operands)
    else
       allocate(C%subs(2))
       C%subs(1)%ptr => left
       C%subs(2)%ptr => right
       
       C%m_dim = left%m_dim
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
    print*, "B%scalar=", B%scalar
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

    print*, "hello!!!"
    A%data  => B%data
    A%node_type = B%node_type
    A%m_dim = B%m_dim
    A%m_shape = B%m_shape
    A%alpha = B%alpha
    A%beta = B%beta

    call copy_node_array(A%subs, B%subs)
    call copy_node_array(A%operands, B%operands)

    A%is_implicit = B%is_implicit
    
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

    !call node_combine(C, res)
    
    call write_graph(C)

    ! call node_destroy(C, ierr)
    
    !write(*, "(I3)") size(C%operands)
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
    
    if(A%node_type == type_data) return

    do i = 1,size(A%subs)
       !combine the left node
       call node_combine(A%subs(i)%ptr, subs)

       !if the node_type does not match its leaf node_type,
       ! not accept the returned list from leaf
       if(A%node_type /= A%subs(i)%ptr%node_type) then
          if(allocated(subs))  deallocate(subs)
          allocate(subs(1))
          subs(1)%ptr => A%subs(i)%ptr
       end if
       
       if(allocated(subs)) then
          cnt = cnt + size(subs)

          call reallocate(res, cnt)
          res(cnt-size(subs)+1 : cnt) = subs
          
          ! if(allocated(res)) then
          !    allocate(tmp(size(res)))
          !    tmp = res
          !    deallocate(res)

          !    allocate(res(cnt))
          !    res(1:size(tmp)) = res
          !    res(size(tmp)+1:cnt) = subs
          !    deallocate(tmp)
          ! else
          !    allocate(res(cnt))
          !    res = subs
          ! endif
       endif

    enddo

    call copy_node_array(A%operands, res)
    
    ! allocate(A%operands(size(res)))
    ! A%operands = res

    ! !combine the left node
    ! call node_combine(A%left, left)

    ! !if the node_type does not match its leaf node_type,
    ! ! not accept the returned list from leaf
    ! if(A%node_type /= A%left%node_type) then
    !    if(allocated(left))  deallocate(left)
    !    allocate(left(1))
    !    left(1)%ptr => A%left
    ! end if

    ! !combine the right node
    ! call node_combine(A%right, right)

    ! if(A%node_type /= A%right%node_type) then
    !    if(allocated(right)) deallocate(right)
    !    allocate(right(1))
    !    right(1)%ptr => A%right
    ! end if

    !combine the result from left and right leaves
    ! allocate(res(size(left) + size(right)))
    ! res(1:size(left)) = left
    ! res(size(left)+1:size(left)+size(right)) = right

    ! allocate(A%operands(size(left) + size(right)))
    ! A%operands = res

    ! write(*,*) "**********************"
    ! if(associated(A%operands)) then
    !    do i = 1, size(A%operands)
    !       write(*, "(Z16.16, A, Z16.16, A, Z16.16)"), &
    !            loc(A%operands(i)%ptr), &
    !            " left=",  loc(A%operands(i)%ptr%left), &
    !            " right=", loc(A%operands(i)%ptr%right)

    !    enddo
    ! end if

    ! deallocate(left, right)
  end subroutine

  recursive subroutine node_destroy(A, ierr)
    implicit none
    type(node), intent(inout), target :: A
    integer :: ierr
    integer i

    if(A%node_type == type_data) return
    
    ! do i = 1, size(A%operands)
    !    write(*, "(Z16.16, A, I4, A, I4)"), &
    !         loc(A%operands(i)%ptr), " = ", A%operands(i)%ptr%node_type, &
    !         " = ", size(A%operands(i)%ptr%operands)
    ! enddo
    
    if(allocated(A%subs)) then
       do i = 1, size(A%subs)
          call node_destroy(A%subs(i)%ptr, ierr)
       enddo
    endif
    
    !destroy the operands
    deallocate(A%operands)
    !A%operands => null()

    !destroy data if the data is implicit
    if(associated(A%data)) then
       if(A%data%is_implicit) then
          call tensor_destroy(A%data, ierr)
          A%data => null()
       end if
    end if
  end subroutine
 
  recursive subroutine write_graph(A, is_root)
    implicit none
    type(node), intent(inout) :: A
    integer :: i
    character(len=100) :: label1, label2
    integer id
    logical,intent(in), optional :: is_root
    logical :: am_i_root
    integer, parameter :: out_unit = 20
    character(len=10) :: op_name
    
    character(len=*), parameter :: &
         format_label = "(I0.1,A,A,A,F4.1,A,F4.1,A)"
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
       open(unit=out_unit,file="graph.dot", &
            action="write", status="replace")
       write(out_unit, *) "digraph G {"
    end if

    select case (A%node_type)
    case (type_data)
       write(op_name, *) "data"
    case (type_plus)
       write(op_name, *) "(+)"
    case (type_minus)
       write(op_name, *) "(-)"
    case (type_mult)
       write(op_name, *) "(*)"
    case (type_divd)
       write(op_name, *) "(/)"
    case default
       write(op_name, *) "func"
    end select
    
    write(out_unit, format_label) id,'[label="',trim(op_name), '\n(',&
         A%alpha, ',',A%beta, ')"];'
    !print*, label1
    if(A%node_type == type_data) then
       return
    end if

    !if(is_root)  print*, "id=", id, am_i_root    
    ! write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"    
    ! write(*, "(A, Z16.16, A, Z16.16, A, Z16.16, I4.3, I3)") &
    !      "A=", loc(A), &
    !      " left=", loc(A%left), " right=", &
    !      loc(A%right), A%node_type, size(A%operands)

    print*, "subs=", size(A%subs)
    write(*, *) "type=", A%node_type

    do i = 1, size(A%subs)
       write(*, "(I0.1, A, Z16.16)"), i," = ", loc(A%subs(i)%ptr)           
       call write_graph(A%subs(i)%ptr,  .false.)
       write(out_unit, format_map), A%id, "->", A%subs(i)%ptr%id, ";"
    end do
    
    if(am_i_root) then
       write(out_unit, *) "}"
       close (out_unit)
    end if

  end subroutine

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

