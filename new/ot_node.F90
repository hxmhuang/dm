#include "type.h"
#:include "type_def.fypp"

module ot_node
  use ot_ref
  use ot_type  
  use ot_tensor
  use ot_vector
  
  interface node_new
     module procedure node_new_from_tensor
     module procedure node_new_from_node
     module procedure node_new_op_binary
     module procedure node_new_op_tensor
     module procedure node_new_op_node     
#:for t in ['integer', 'real', 'real8']
     module procedure node_new_from_${t}$
#:endfor
  end interface node_new
  
  !private :: shape 
  interface shape
     module procedure node_shape
  end interface shape

  ! interface dim
  !    module procedure node_dim
  ! end interface dim

  interface box
     module procedure node_box
  end interface box

  interface lbox
     module procedure node_lbox
  end interface lbox

  interface reallocate
     module procedure reallocate1,reallocate2
  end interface reallocate
  
  character(len=8), dimension(${L[len(L)-1][1]}$) :: op_names

  interface disp_info
     module procedure disp_info_node     
  end interface disp_info

  interface disp
     module procedure disp_node
  end interface disp
  
contains

  subroutine init_node(ierr)
    implicit none
    integer, intent(out) :: ierr
    
#:for e in L
    op_names(${e[1]}$) = "${e[3]}$"
#:endfor

  end subroutine
  
  !***********************
  !> node constructers
  !***********************
  subroutine node_new_from_tensor(B, A) 
    implicit none
    type(tensor), intent(in), target :: A
    type(node),   intent(out) :: B

    !call tensor_ensure_valid(A)
    TENSOR_ENSURE_VALID(A)
    
    B%id = get_global_id()
    B%m_shape = A%m_shape
    call assign_ptr(B%data, A)
    B%node_type = type_data
  end subroutine

  !copy node A to node B
  subroutine node_new_from_node(B, A) 
    implicit none
    type(node), intent(in)   :: A
    type(node), intent(out)  :: B
    integer i

    B%id = A%id
    B%m_shape = A%m_shape
    B%alpha = A%alpha
    B%beta = A%beta

    if(associated(A%data)) &
         call assign_ptr(B%data, A%data)
    
    B%scalar = A%scalar
    B%args = A%args
    B%node_type = A%node_type
    B%ref = A%ref
    call copy_node_ptr(B%operands, A%operands)
  end subroutine

  subroutine node_new_op_binary(C, op_type, left, right)
    implicit none
    integer, intent(in) :: op_type
    type(node), intent(in), target :: left, right
    type(node), intent(out) :: C
    real(8) :: x

    if((.not. is_scalar(left)) .and. &
         (.not. is_scalar(right))) then
       
       call assert(.not. any(left%m_shape /= right%m_shape), &
            __FILE__, __LINE__, &
            "shape of the left and right leaf does not match.")
    endif
    
    if(is_scalar(left)) then
       x = left % scalar
       select case(op_type)
       case (type_plus) ! x + A
          call node_new(C, right)
          C%alpha = x + C%alpha
       case (type_minus) ! x - A
          call node_new(C, right)
          C%alpha = x - C%alpha
          C%beta  = -C%beta
       case (type_mult) ! x .* A
          call node_new(C, right)
          C%alpha = x * C%alpha
          C%beta  = x * C%beta
       case (type_divd) ! x ./ A
          C%m_shape = right%m_shape
          C%beta = x
          C%node_type = type_rcp
          ! allocate(C%operands(1))
          ! C%operands(1)%ptr => right
          call push_back(C%operands, right)
       end select
    else if(is_scalar(right)) then
       x = right % scalar
       call node_new(C, left)      
       select case(op_type)
       case (type_plus) ! A + x
          C%alpha = C%alpha + x
          C%beta  = C%beta
       case (type_minus) ! A - a
          C%alpha = C%alpha - x
          C%beta  = C%beta
       case (type_mult) ! A .* a
          C%alpha = C%alpha * x
          C%beta  = C%beta  * x
       case (type_divd)
          C%alpha = C%alpha / x
          C%alpha = C%beta  / x
       end select
    else
       allocate(C%operands(2))
       ! C%operands(1)%ptr => left
       ! C%operands(2)%ptr => right
       ! write(*, "(A, Z16.16)"), "left=", loc(left)
       ! write(*, "(A, Z16.16)"), "right=", loc(right)
       call assign_ptr(C%operands(1)%ptr, left)
       call assign_ptr(C%operands(2)%ptr, right)
       if(is_scalar(left)) then
          C%m_shape = right%m_shape
       else
          C%m_shape = left%m_shape
       endif
       C%node_type = op_type
    end if
  end subroutine

#:for t in ['integer', 'real', 'real8']
  !> create a node from a ${t}$ scalar
  subroutine node_new_from_${t}$(B, scalar) 
    implicit none
    type(${t}$), intent(in), target :: scalar
    type(node), intent(out) :: B

    B%id = get_global_id()
    B%m_shape = (/1, 1, 1/)
    B%scalar = real(scalar, 8)
    B%node_type = type_scalar
  end subroutine
#:endfor

  !> operation on a tensor
  subroutine node_new_op_tensor(dst, op_type, src)
    implicit none
    type(node), intent(out) :: dst
    integer :: op_type
    type(tensor), intent(in) :: src
    type(node) :: dst1

    !call tensor_ensure_valid(src)
    TENSOR_ENSURE_VALID(src)

    call node_new(dst1, src)
    
    call node_new(dst, op_type, dst1)

    !print*, trim(op_names(op_type))
  end subroutine

  !> operation on a node
  subroutine node_new_op_node(dst, op_type, src)
    implicit none
    type(node), intent(out) :: dst
    integer :: op_type
    type(node), intent(in) :: src
    type(node), pointer :: p
    
    dst%node_type = op_type
    dst%m_shape   = src%m_shape
    !dst%local_block = src%local_block

    allocate(p)
    call node_new(p, src)

    call push_back(dst%operands, p)    
    ! allocate(dst%operands(1))
    ! dst%operands(1)%ptr => p
  end subroutine

  subroutine disp_info_node(o, prefix)
    implicit none
    type(node), intent(in) :: o
    character(len=*), optional, intent(in) :: prefix
    integer :: i
    
    write(*,*) ""
    if(present(prefix)) then
       write(*, "(A)") prefix
    else
       write(*, "(A)") "ans = "
    endif

    write(*, "(4X, A, I0.1)") "node id : ", o%id
#ifdef DEBUG
    write(*, "(4X, A, Z16.16)") "obj addr : 0x", loc(o)
#endif
    write(*, "(4X, A, A)") "node type : ", op_names(o%node_type)    
    write(*, "(4X, A)", advance="no") "shape : ["
    do i = 1, size(o%m_shape)
       write(*, "(I0.1)", advance="no") o%m_shape(i)
       if(i < size(o%m_shape)) write(*, "(A)", advance="no") "x"
    enddo
    write(*, "(A)") "]"

    write(*, "(4X, A, F8.4)") "alpha : ", o%alpha
    write(*, "(4X, A, F8.4)") "beta  : ", o%beta
    if(is_scalar(o)) write(*, "(A, F0.8)") "scalar  : ", o%scalar
    write(*, "(4X, A, 10F8.4)") "args : ", o%args(1:5)

    if(.not. is_data(o)) then
       write(*, "(4X, A)") "operands : "
       do i = 1, size(o%operands)
          write(*, "(6X, A, Z16.16)") "0x", loc(o%operands(i)%ptr)
       enddo
    else
       write(*, "(4X, A, Z16.16)") "data : 0x", loc(o%data)
    endif

    write(*, "(A)") ""
  end subroutine

  subroutine disp_node(A, msg)
    implicit none
    type(node) :: A
    character(len=*),intent(in),optional :: msg
    if(present(msg)) then
       call disp(A%data, msg)
    else
       call disp(A%data)
    end if
  end subroutine

  
  function node_dim(o) result(res)
    implicit none
    type(node), intent(in) :: o
    integer :: res
    
  end function

  function node_shape(o) result(res)
    implicit none
    type(node), intent(in) :: o
    integer :: res(3)
    res = o%m_shape
  end function

  function node_lbox(o) result(res)
    implicit none
    type(node), intent(in) :: o
    type(box_info) :: res

    if(associated(o%data)) then
       res = lbox(o%data)
    else
       !this is an invalid box, means unknow
       res%starts = 0
       res%ends   = -1
    endif

  end function

  function node_box(o) result(res)
    implicit none
    type(node), intent(in) :: o
    type(box_info) :: res

    res%starts = 0
    res%ends = res%starts + shape(o) - 1

  end function

  function is_scalar(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    res = (A%node_type == type_scalar)
  end function

  function is_data(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    res = (A%node_type == type_scalar) .or. &
         (A%node_type == type_data) .or. &
         (A%node_type == type_ref) 
  end function

  function is_ref(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    res = (A%node_type == type_ref) 
  end function
  
  function is_arithmetic(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    res = (A%node_type == type_plus) .or. &
         (A%node_type == type_minus) .or. &
         (A%node_type == type_mult)  .or. &
         (A%node_type == type_divd)
  end function

  function need_eval(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    integer :: i

    res = any(A%args /= 0)
    ! do i = 1, size(A%args)
    !    if(A%args(i) .ne. 0.) then
    !       res = .true.
    !       return
    !    endif
    ! enddo
    
    res =res .or. (A%alpha /= 0 .or. &
         A%beta /= 1.0 .or. &
         (.not. associated(A%data)))
    
  end function

  subroutine copy_node_ptr(dst, src)
    implicit none
    type(node_ptr), allocatable, intent(in) :: src(:)
    type(node_ptr), allocatable, intent(inout) :: dst(:)

    if(allocated(src)) then
       if(allocated(dst)) deallocate(dst)
       allocate(dst(size(src)))
       dst = src
    end if
  end subroutine

  subroutine node_add_operand(t, o)
    implicit none
    type(node), intent(inout) :: t
    type(node), intent(in), pointer :: o

    call push_back(t%operands, o)
  end subroutine

  !> reallocate an array while keeping its data
  subroutine reallocate1(arr, len)
    implicit none
    type(node_ptr), allocatable, intent(inout) :: arr(:)
    type(node_ptr), allocatable :: tmp(:)
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
    type(node_ptr), pointer, intent(inout) :: arr(:)
    type(node_ptr), pointer :: tmp(:)
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

  recursive subroutine node_destroy(A, ierr)
    implicit none
    type(node), intent(inout), target :: A
    integer :: ierr
    integer i

    ! A%ref_cnt = A%ref_cnt - 1
    ! if(A%ref_cnt < 0) A%ref_cnt = 0
    ! if(A%ref_cnt <= 0) return
    
    if(A%node_type == type_data .or. &
         A%node_type == type_scalar) return
    
    if(allocated(A%operands)) then
       do i = 1, size(A%operands)
          call node_destroy(A%operands(i)%ptr, ierr)
       enddo
    endif
    
  end subroutine

  recursive subroutine disp_tree(o)
    implicit none
    type(node), intent(in) :: o
    integer :: i
    
    if(is_data(o)) return
    call disp_info(o)
    
    do i = 1, size(o%operands)
       call disp_tree(o%operands(i)%ptr)
    end do
    
  end subroutine
end module
