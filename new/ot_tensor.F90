#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 
#:include "type_def.fypp"
#include "type.h"

module ot_tensor
  use ot_common
  use ot_data
  use ot_type  
  use ot_ref

  !here we use a trick to extract range information from an empty array
  !see find_range function
  !integer, allocatable :: r(:)
  interface
     subroutine reg_func(p, id) &
          bind(C, name="reg_func")
       implicit none
       C_POINTER :: p
       integer :: id
     end subroutine
  end interface

  interface shape
     module procedure tensor_shape     
  end interface shape

  ! interface dim
  !    module procedure tensor_dim
  ! end interface dim

  interface box
     module procedure tensor_box
  end interface box

  interface lbox
     module procedure tensor_lbox
  end interface lbox

  interface disp_info
     module procedure disp_info_tensor
  end interface disp_info

  interface tensor_new
     module procedure tensor_new_from_data
  end interface tensor_new

  interface disp
     module procedure disp_tensor
  end interface disp

  interface destroy
     module procedure tensor_destroy
  end interface destroy

  interface is_valid
     module procedure is_valid_tensor
  end interface is_valid

  interface assign_ptr
     module procedure assign_tensor_ptr
  end interface assign_ptr

  interface release_ptr
     module procedure release_tensor_ptr
  end interface release_ptr

#:set t = 'tensor'
#:include "ot_vector.if"
#:del t
  
contains

#:set t = 'tensor'
#:include "ot_vector.inc"
#:del t
  
  subroutine init_tensor(ierr)
    implicit none
    integer, intent(out) :: ierr

    !initialize data module
    call init_data(ierr)

  end subroutine

  
  subroutine assign_tensor_ptr(p, o)
    implicit none
    type(tensor), pointer,intent(out) :: p
    type(tensor), target, intent(in)  :: o
    integer :: cnt
    p => o
    ! add reference counter
    cnt = inc_ref_cnt(p)
  end subroutine
  
  !>delete tensor pointer
  subroutine release_tensor_ptr(p)
    implicit none
    type(tensor), pointer,intent(inout) :: p
    integer :: ierr
    
    if(dec_ref_cnt(p) <= 0) then
       if(is_rvalue(p)) then
          call destroy(p, ierr)
          deallocate(p)
          p => null()
       endif
    endif

  end subroutine

  !>destroy the tensor
  subroutine tensor_destroy(A, ierr)
    implicit none
    type(tensor), intent(inout) :: A
    integer, intent(out) :: ierr
    integer :: cnt

    call data_destroy(A%data, ierr)
  end subroutine 

  subroutine tensor_new_from_data(t, data)
    implicit none
    type(tensor), intent(out) :: t
    Vec, intent(in) :: data

    call assert(data /= 0, __FILE__, __LINE__, &
         "data can not be null.")
    
    call tensor_bind_data(t, data)
  end subroutine

  subroutine tensor_bind_data(o, data)
    implicit none
    type(tensor), intent(inout) :: o
    Vec, intent(in) :: data
    integer :: d_shape(3)
    integer :: dx, dy, dz

    o%data = data    
    call petsc_get_corners(data, o%local_block)
    call petsc_data_info(data, dx, dy, dz)
    o%m_shape = (/dx, dy, dz/)
  end subroutine
  
  function find_range(ra) result(res)
    implicit none
    integer, intent(inout) :: ra(:)
    integer :: res(2)

    !print*, size(range)
    res(1) = (loc(ra) - loc(r)) / 4
    res(2) = ((loc(ra) + size(ra) * 4) - loc(r) - 4)/4

  end function
  
  function tensor_shape(o) result(res)
    implicit none
    type(tensor), intent(in) :: o
    integer :: res(3)
    res = o%m_shape
  end function
  
  function tensor_box(o) result(res)
    implicit none
    type(tensor), intent(in) :: o
    type(box_info) :: res
    integer :: s(3)

    s = shape(o)
    
    res%rx = r(0, s(1))
    res%ry = r(0, s(2))
    res%rz = r(0, s(3))    
  end function

  function tensor_lbox(o) result(res)
    implicit none
    type(tensor), intent(in) :: o
    type(box_info) :: res
    integer :: xs, ys, zs, xl, yl, zl

    call petsc_get_corners(o%data, res)
  end function

  subroutine disp_info_tensor(o, prefix)
    implicit none
    type(tensor), intent(in) :: o
    character(len=*), optional, intent(in) :: prefix
    integer :: i
    integer :: dim
    
    if(get_rank() == 0) then
       write(*,*) ""
       if(present(prefix)) then
          write(*, "(A)") prefix
       else
          write(*, "(A)") "ans = "
       endif
#ifdef DEBUG
       write(*, "(4X, A, Z16.16)") "obj addr : 0x", loc(o)
       write(*, "(4X, A, I3.1)") "ref count :", o%ref_cnt
#endif
       call disp_shape(o%m_shape, '4X')
       ! dim = find_dim(o%m_shape)
       ! write(*, "(4X, A)", advance="no") "shape : ["
       ! if(dim == 0) then
       !    write(*, "(I0.1)", advance="no") 0
       ! else
       !    do i = 1, dim
       !       write(*, "(I0.1)", advance="no") o%m_shape(i)
       !       if(i < dim) &
       !            write(*, "(A)", advance="no") "x"
       !    enddo
       ! endif
       write(*, "(A)") "]"
       
       write(*, "(4X, A, L1)") "is_field    : ", o%is_field
       if(o%is_field) &
            write(*, "(4X, A, I0.2)") &
            "grid_pos : ", o%grid_pos
       write(*, "(4X, A, Z16.16)") "data : 0x", o%data
    end if
  end subroutine


  function is_valid_tensor(t) result(res)
    implicit none
    type(tensor) :: t
    logical :: res
    
    if(any(shape(t) <= 0) &
         .or. t%data /=0) then
       res = .true.
    else
       res = .false.
    end if
  end function
  
  subroutine disp_tensor(objA, prefix)
#include "petsc.h"            
    type(tensor), intent(in),target :: objA
    type(tensor), pointer :: A
    character(len=*), optional, intent(in) :: prefix
    real(kind=8), pointer :: x1(:), x2(:,:), x3(:,:,:)
    integer :: i, ierr, rank
    integer :: dim
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    A => objA
    !call VecView(objA%data, PETSC_VIEWER_STDOUT_WORLD,ierr)

    if(rank == 0) then        
       write(*,*) ""
       if(present(prefix)) then
          write(*, "(A)") prefix
       else
          write(*, "(A)") "ans = "
       endif

       call disp_shape(A%m_shape, '4X')
       ! dim = find_dim(A%m_shape)
       ! write(*, "(4X, A)", advance="no") "shape : ["
       ! do i = 1, dim
       !    write(*, "(I0.1)", advance="no") A%m_shape(i)
       !    if(i < dim) write(*, "(A)", advance="no") "x"
       ! enddo
       ! write(*, "(A)") "]"
    endif
    
    if(A%data /= 0) then
       call petsc_print(A%data)
    end if
  end subroutine

  subroutine tensor_set_shape(A, m_shape)
    implicit none
    type(tensor), intent(out) :: A    
    integer, intent(in) :: m_shape(:)
    
    select case (size(m_shape))
    case (1)
       A%m_shape(1) = m_shape(1)
       A%m_shape(2) = 1
       A%m_shape(3) = 1       
    case (2)
       A%m_shape(1) = m_shape(1)
       A%m_shape(2) = m_shape(2)
       A%m_shape(3) = 1       
    case (3)
       A%m_shape(1) = m_shape(1)
       A%m_shape(2) = m_shape(2)
       A%m_shape(3) = m_shape(3)
    end select
    
  end subroutine


#:for op in L
  #:if op[1] >= 50
subroutine ${op[2]}$_tensors(res, tensor_operands, &
     ops_alpha, ops_beta, args) 
    implicit none
#include "petsc.h"
    type(tensor), intent(inout) :: res
    type(tensor_ptr), intent(in), allocatable  :: tensor_operands(:)
    Vec,allocatable :: vec_operands(:)
    real(8), intent(in) :: ops_alpha(:), ops_beta(:)!, alpha, beta
    real(8), intent(in) :: args(10)
    integer :: n
    integer :: ierr, i

    call assert(allocated(tensor_operands), &
         __FILE__, __LINE__, &
         "tensor_operands must be allocated!")

    n = size(tensor_operands)    
    allocate(vec_operands(n))

    !print*, "loc(res)=", loc(res)
    !print*, "loc(ptr)=", loc(tensor_operands(1)%ptr)    
    !call VecView(tensor_operands(1)%ptr%data, PETSC_VIEWER_STDOUT_WORLD,ierr)
    !write(*, "(A, Z16.16)"), "ptr_data=", tensor_operands(1)%ptr%data
    !call test_tensor(res, tensor_operands(1)%ptr%data)    
    !write(*, "(A, Z16.16)"), "ptr_data=", tensor_operands(1)%ptr%data
    
    do i = 1, n
       vec_operands(i) = tensor_operands(i)%ptr%data
    enddo

#ifdef DEBUG
    print*, "calling data_${op[2]}$"
#endif
    
    ! call VecView(res%data, PETSC_VIEWER_STDOUT_WORLD,ierr)
    ! call VecView(vec_operands(1), PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    call data_${op[2]}$(res%data, vec_operands, &
         ops_alpha, ops_beta, args, ierr)

    ! call VecView(res%data, PETSC_VIEWER_STDOUT_WORLD,ierr)
    
  end subroutine
  #:endif
#:endfor

  subroutine tensor_copy_structure(dst, src)
    implicit none
    type(tensor), intent(inout) :: dst
    type(tensor), intent(in) :: src

    dst%m_shape  = src%m_shape
    dst%is_field = src%is_field
    dst%grid_pos = src%grid_pos
  end subroutine
  
  !> naively copy the data member and pointers from src to dst
  subroutine tensor_copy(dst, src)
    implicit none
    type(tensor), intent(inout) :: dst
    type(tensor), intent(in)    :: src
    integer :: ierr
    
    call tensor_copy_structure(dst, src)
    call data_destroy(dst%data,  ierr)    
    dst%data = src%data
  end subroutine

  !> deep copy tensor, all data are copied from dst to src
  subroutine tensor_deep_copy(dst, src)
    implicit none
    type(tensor), intent(out) :: dst
    type(tensor), intent(in)  :: src
    integer :: ierr

    call tensor_copy_structure(dst, src)
    call data_destroy(dst%data,  ierr) !safe destroy
    call data_clone(dst%data, src%data, ierr)
  end subroutine 

  !> duplicate tensor structure, the data is allocated but not copied
  subroutine tensor_duplicate(dst, src)
    implicit none
    type(tensor), intent(out) :: dst
    type(tensor), intent(in), pointer  :: src
    integer :: ierr

    !write(*, "(A, Z16.16)"), "src_data=", src%data
    
    call tensor_copy_structure(dst, src)
    call data_destroy(dst%data,  ierr) !safe destroy
    call data_duplicate(dst%data,  src%data,  ierr)
  end subroutine 

  subroutine test_tensor(dst, v)
    implicit none
#include "petsc.h"
    type(tensor), intent(out) :: dst
    !type(tensor), pointer  :: src
    Vec, intent(in) :: v
    integer :: ierr

    !write(*, "(A, Z16.16)"), "src_data=", src%data
    write(*, "(A, Z16.16)"), "v=", v
    
    ! call tensor_copy_structure(dst, src)
    ! call data_destroy(dst%data,  ierr) !safe destroy
    ! call data_duplicate(dst%data,  src%data,  ierr)
  end subroutine

  subroutine slice_tensors(dst, src, ref)
    implicit none
#include "petsc.h"
    type(tensor), intent(inout) :: dst
    type(tensor), intent(in)    :: src
    type(ref_info), intent(in)  :: ref
    Vec :: dst_data
    integer :: ierr
    
    !call data_get_sub(dst%data, src%data, box)
    call data_get_sub(dst_data, src%data, ref)
    call tensor_new(dst, dst_data)

    !call VecView(dst%data, PETSC_VIEWER_STDOUT_WORLD, ierr)
  end subroutine
end module 

