#include "type.h"
#:include "type_def.fypp"

module ot_node
  use ot_ref
  use ot_type  
  use ot_tensor
  
  interface node_new
     module procedure node_new_from_tensor
     module procedure node_new_from_node
     module procedure node_new_op_binary
     !module procedure node_new_op_tensor
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

  interface disp_info
     module procedure disp_info_node     
  end interface disp_info

  interface destroy
     module procedure node_destroy
  end interface destroy
  
  interface is_valid
     module procedure is_valid_node
  end interface is_valid

  interface assign_ptr
     module procedure assign_node_ptr
  end interface assign_ptr

  interface release_ptr
     module procedure release_node_ptr
  end interface release_ptr

#:set t = 'node'
#:include "ot_vector.if"
#:del t
  
  character(len=8), dimension(${L[len(L)-1][1]}$) :: op_names

  type op_desc
     character(len=8) :: op_name
     character(len=1) :: ntype ! unary/binary
     character(len=1) :: etype ! element-wise type
     character(len=20) :: expr_form !expression form
  end type op_desc

  type(op_desc), dimension(${L[len(L)-1][1]}$) :: op_desc_list
contains

#:set t = 'node'
#:include "ot_vector.inc"
#:del t
  
  !> initialize the node module
  subroutine init_node(ierr)
    implicit none
    integer, intent(out) :: ierr
    integer :: i
    
#:for e in L
    op_names(${e[1]}$) = "${e[3]}$"
    
    op_desc_list(${e[1]}$)%op_name   = "${e[3]}$"
    op_desc_list(${e[1]}$)%ntype     = "${e[4]}$"
    op_desc_list(${e[1]}$)%etype     = "${e[5]}$"
    op_desc_list(${e[1]}$)%expr_form = "${e[6]}$"        
#:endfor

    ! do i = 1, size(op_desc_list)
    !    print*, 'i=', i, 'etype=', op_desc_list(i)%etype
    ! enddo
  end subroutine

  !> is element-wise
  function is_ew(o) result(res)
    type(node), intent(in) :: o
    logical :: res

    res = (op_desc_list(o%node_type)%etype == 'T')
  end function
  
  subroutine assign_node_ptr(p, o)
    implicit none
    type(node), pointer,intent(out) :: p
    type(node), target, intent(in)  :: o
    integer :: cnt
    p => o
    ! add reference counter
    cnt = inc_ref_cnt(p)
  end subroutine

  !> delete node pointer
  recursive subroutine release_node_ptr(p)
    implicit none
    type(node), pointer,intent(inout) :: p
    integer :: ierr
    
    if(dec_ref_cnt(p) <= 0) then
       if(is_rvalue(p)) then
          ! call disp_info(p, 'deleting node = ')
          ! if(associated(p%data)) then
          !    print*, 'var_type = ', p%data%var_type
          ! end if
          call destroy(p, ierr)
          deallocate(p)
       end if
    endif
  end subroutine

  !>destroy node A
  recursive subroutine node_destroy(A, ierr)
    implicit none
    type(node), intent(inout) :: A
    integer :: ierr
    integer i, cnt

    if(associated(A%data)) then
       call release_ptr(A%data)
    end if

    if(allocated(A%operands)) then
       do i = 1, size(A%operands)
          call release_ptr(A%operands(i)%ptr)
       enddo
       deallocate(A%operands)
    endif
  end subroutine
  
  !> node constructers
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

  !>copy node A to node B
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

  !> create a node from binary operation
  subroutine node_new_op_binary(C, op_type, left, right)
    implicit none
    integer, intent(in) :: op_type
    type(node), intent(in), pointer :: left, right
    type(node), intent(out), pointer :: C
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
          call assign_ptr(C, right)
          C%alpha = x + C%alpha
       case (type_minus) ! x - A
          call assign_ptr(C, right)          
          C%alpha = x - C%alpha
          C%beta  = -C%beta
       case (type_mult) ! x .* A
          call assign_ptr(C, right)          
          C%alpha = x * C%alpha
          C%beta  = x * C%beta
       case (type_divd) ! x ./ A
          allocate(C)
          C%id = get_global_id()
          C%m_shape = right%m_shape
          C%beta = x
          C%node_type = type_rcp
          call push_back(C%operands, right)
       end select
    else if(is_scalar(right)) then
       x = right % scalar
       call assign_ptr(C, left)
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
       allocate(C)
       allocate(C%operands(2))
       
       call assign_ptr(C%operands(1)%ptr, left)
       call assign_ptr(C%operands(2)%ptr, right)
       
       if(is_scalar(left)) then
          C%m_shape = right%m_shape
       else
          C%m_shape = left%m_shape
       endif
       
       C%node_type = op_type
       
       !assign new node id
       C%id = get_global_id()
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

  !> operation on a node
  subroutine node_new_op_node(dst, op_type, src)
    implicit none
    type(node), intent(out) :: dst
    integer :: op_type
    type(node), intent(in), target :: src
    type(node), pointer :: p
    
    dst%node_type = op_type
    dst%m_shape   = src%m_shape

    p => src
    ! allocate(p)
    ! call node_new(p, src)

    call push_back(dst%operands, p)
    
    dst%id = get_global_id()
  end subroutine

  subroutine disp_info_node(o, prefix)
    implicit none
    type(node), intent(in) :: o
    character(len=*), optional, intent(in) :: prefix
    integer :: i, dim
    
    write(*,*) ""
    if(present(prefix)) then
       write(*, "(A)") prefix
    else
       write(*, "(A)") "ans = "
    endif

    write(*, "(4X, A, I0.1)") "node id : ", o%id
#ifdef DEBUG
    write(*, "(4X, A, Z16.16)") "obj addr : 0x", loc(o)
    write(*, "(4X, A, I3.1)") "ref count :", o%ref_cnt
    if(associated(o%data)) then
       write(*, "(4X, A, Z16.16, A, I3.1)") "data addr :", loc(o%data), &
            ", ref_cnt = ", o%data%ref_cnt
    end if
#endif
    write(*, "(4X, A, A)") "node type : ", op_names(o%node_type)
    call disp_shape(o%m_shape, '4X')
    !write(*, "(4X, A)", advance="no") "shape : ["
    ! dim = find_dim(o%m_shape)
    ! if(dim == 0) then
    !    write(*, "(I0.1)", advance="no") 0
    ! else
    !    do i = 1, dim
    !       write(*, "(I0.1)", advance="no") o%m_shape(i)
    !       if(i < dim) write(*, "(A)", advance="no") "x"
    !    enddo
    ! endif
    ! write(*, "(A)") "]"

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
       res%rx = r(0, -1)
       res%ry = r(0, -1)
       res%rz = r(0, -1)       
    endif

  end function

  function node_box(o) result(res)
    implicit none
    type(node), intent(in) :: o
    type(box_info) :: res

    res%rx = r(0, -1)
    res%ry = r(0, -1)
    res%rz = r(0, -1)       

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

  function is_valid_node(t) result(res)
    implicit none
    type(node) :: t
    logical :: res
    
    if(any(shape(t) <= 0)) then
       res = .true.
    else
       res = .false.
    end if
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

  recursive subroutine disp_tree(o)
    implicit none
    type(node), intent(in) :: o
    integer :: i

    if(get_rank() == 0) then
       call disp_info(o)

       if(is_data(o)) return
       
       do i = 1, size(o%operands)
          call disp_tree(o%operands(i)%ptr)
       end do
    endif
  end subroutine
end module
