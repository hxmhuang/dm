#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 
#:include "type_def.fypp"
#include "type.h"

module ot_array
  use ot_common
  use ot_type  
  use ot_ref
  use ot_petsc
  
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

  interface array_new
     module procedure array_new_from_data
  end interface array_new

  interface disp
     module procedure disp_array
  end interface disp

  interface shape
     module procedure array_shape     
  end interface shape

  interface box
     module procedure array_box
  end interface box

  interface lbox
     module procedure array_lbox
  end interface lbox

  interface is_scalar
     module procedure is_scalar_array
  end interface is_scalar

  interface is_valid
     module procedure is_valid_array
  end interface is_valid

  interface disp_info
     module procedure disp_info_array
  end interface disp_info

  
  interface array_seq
     module procedure array_new_seq_with_shape1
     module procedure array_new_seq_with_shape2
     module procedure array_new_seq_with_shape3
     module procedure array_new_seq_with_data
     module procedure array_new_scalar     
  end interface array_seq
  
  interface destroy
     module procedure array_destroy
  end interface destroy

  interface assign_ptr
     module procedure assign_array_ptr
  end interface assign_ptr

  interface bind_ptr
     module procedure bind_array_ptr
  end interface bind_ptr
  
  interface release_ptr
     module procedure release_array_ptr
  end interface release_ptr

  
#:set t = 'array'
#:include "ot_vector.if"
#:del t
  
contains

#:set t = 'array'
#:include "ot_vector.inc"
#:del t
  
  subroutine init_array(ierr)
    implicit none
    integer, intent(out) :: ierr

    !initialize data module
    !call init_data(ierr)

  end subroutine
  
  subroutine assign_array_ptr(p, o)
    implicit none
    type(array), pointer,intent(out) :: p
    type(array), target, intent(in)  :: o
    integer :: cnt
    p => o
    ! add reference counter
    cnt = inc_ref_cnt(p)
  end subroutine

  subroutine bind_array_ptr(p, o)
    implicit none
    type(array), pointer,intent(out) :: p
    type(array), target, intent(in)  :: o
    p => o
  end subroutine
  
  !>delete array pointer
  subroutine release_array_ptr(p)
    implicit none
    type(array), pointer,intent(inout) :: p
    integer :: ierr

    if(.not. associated(p)) return
    
    if(dec_ref_cnt(p) <= 0) then
       if(is_rvalue(p)) then
          call destroy(p, ierr)
          deallocate(p)
          p => null()
       endif
    endif

  end subroutine

  !>destroy the array
  subroutine array_destroy(A, ierr)
    implicit none
    type(array), intent(inout) :: A
    integer, intent(out) :: ierr
    integer :: cnt

    call petsc_data_destroy(A%data, ierr)
  end subroutine 

  subroutine array_new_from_data(t, data)
    implicit none
    type(array), intent(out) :: t
    Vec, intent(in) :: data

    call assert(data /= 0, __FILE__, __LINE__, &
         "data can not be null.")
    
    call array_bind_data(t, data)
  end subroutine

  !> overload function of array_seq
  subroutine array_new_seq_with_shape1(t, m, ierr)
    implicit none
    type(array), intent(out) :: t
    integer, intent(in)  :: m
    integer, intent(out) :: ierr
    Vec :: vec_data

    call petsc_create3d_seq_with_shape(&
         vec_data, m, 1, 1, ierr)
    call array_bind_data(t, vec_data)
  end subroutine

  !> overload function of array_seq  
  subroutine array_new_seq_with_shape2(t, m,n, ierr)
    implicit none
    type(array), intent(out) :: t
    integer, intent(in)  :: m, n
    integer, intent(out) :: ierr
    Vec :: vec_data

    call petsc_create3d_seq_with_shape(&
         vec_data, m, n, 1, ierr)
    call array_bind_data(t, vec_data)
  end subroutine

  !> overload function of array_seq  
  subroutine array_new_seq_with_shape3(t, m,n,k, ierr)
    implicit none
    type(array), intent(out) :: t
    integer, intent(in)  :: m, n, k
    integer, intent(out) :: ierr
    Vec :: vec_data

    call petsc_create3d_seq_with_shape(&
         vec_data, m, n, k, ierr)
    call array_bind_data(t, vec_data)
  end subroutine

  !> overload function of array_seq  
  subroutine array_new_seq_with_data(t, data, ierr)
    implicit none
    type(array), intent(out) :: t
    real(8), intent(in)  :: data(:,:,:)
    integer, intent(out) :: ierr
    integer :: s(3)
    Vec :: vec_data

    s = 0 
    s = shape(data)

    where(s == 0) s = 1
    
    call petsc_create3d_seq_with_data(&
         vec_data, s(1), s(2), s(3), data, ierr)
    call array_bind_data(t, vec_data)
  end subroutine

  !> overload function of array_seq  
  subroutine array_new_scalar(t, data, ierr)
    implicit none
    type(array), intent(out) :: t
    real(8), intent(in)  :: data
    real(8) :: data3d(1,1,1)
    integer, intent(out) :: ierr
    integer :: s(3)
    Vec :: vec_data

    data3d(1,1,1) = data
    call petsc_create3d_seq_with_data &
         (vec_data, 1, 1, 1, data3d, ierr)
    call array_bind_data(t, vec_data)
  end subroutine

  subroutine array_bind_data(o, data)
    implicit none
    type(array), intent(inout) :: o
    Vec, intent(in) :: data
    integer :: d_shape(3)
    integer :: dx, dy, dz

    o%data = data    
    call petsc_get_corners(data, o%local_block)
    call disp(o%local_block, 'block = ')
    call petsc_data_info(data, dx, dy, dz)
    o%m_shape = (/dx, dy, dz/)
  end subroutine
  
  ! function find_range(ra) result(res)
  !   implicit none
  !   integer, intent(inout) :: ra(:)
  !   integer :: res(2)

  !   !print*, size(range)
  !   res(1) = (loc(ra) - loc(r)) / 4
  !   res(2) = ((loc(ra) + size(ra) * 4) - loc(r) - 4)/4

  ! end function
  

  
  subroutine disp_array(objA, prefix)
#include "petsc.h"            
    type(array), intent(in),target :: objA
    type(array), pointer :: A
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

#ifdef DEBUG
       call disp(A%local_block, &
            'local_block = ', '4X', '6X')
#endif
    endif
    
    if(A%data /= 0) then
       call petsc_print(A%data)
    end if
  end subroutine

  subroutine array_set_shape(A, m_shape)
    implicit none
    type(array), intent(out) :: A    
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

  subroutine array_copy_structure(dst, src)
    implicit none
    type(array), intent(inout) :: dst
    type(array), intent(in) :: src

    dst%m_shape  = src%m_shape
    dst%is_field = src%is_field
    dst%grid_pos = src%grid_pos
  end subroutine
  
  !> naively copy the data member and pointers from src to dst
  subroutine array_copy(dst, src)
    implicit none
    type(array), intent(inout) :: dst
    type(array), intent(in)    :: src
    integer :: ierr
    
    call array_copy_structure(dst, src)
    call petsc_data_destroy(dst%data,  ierr)    
    call array_bind_data(dst, src%data)
    
    !dst%data = src%data    
  end subroutine

  !> deep copy array, all data are copied from dst to src
  subroutine array_deep_copy(dst, src)
    implicit none
    type(array), intent(inout) :: dst
    type(array), intent(in)  :: src
    integer :: ierr

    call array_copy_structure(dst, src)
    call petsc_data_destroy(dst%data,  ierr) !safe destroy
    call petsc_data_clone(dst%data, src%data, ierr)
  end subroutine 

  !> duplicate array structure, the data is allocated but not copied
  subroutine array_duplicate(dst, src)
    implicit none
    type(array), intent(inout) :: dst
    type(array), intent(in), pointer  :: src
    integer :: ierr

    !write(*, "(A, Z16.16)"), "src_data=", src%data
    
    call array_copy_structure(dst, src)
    call petsc_data_destroy(dst%data,  ierr) !safe destroy
    call petsc_data_duplicate(dst%data,  src%data,  ierr)
  end subroutine 

!   subroutine slice_arrays(dst, src, ref)
!     implicit none
! #include "petsc.h"
!     type(array), intent(inout) :: dst
!     type(array), intent(in)    :: src
!     type(ref_info), intent(in)  :: ref
!     Vec :: dst_data
!     integer :: ierr
    
!     !call data_get_sub(dst%data, src%data, box)
!     call data_get_sub(dst_data, src%data, ref)
!     call array_new(dst, dst_data)

!     !call VecView(dst%data, PETSC_VIEWER_STDOUT_WORLD, ierr)
  !   end subroutine


    function array_shape(o) result(res)
    implicit none
    type(array), intent(in) :: o
    integer :: res(3)
    res = o%m_shape
  end function
  
  function array_box(o) result(res)
    implicit none
    type(array), intent(in) :: o
    type(box_info) :: res
    integer :: s(3)

    s = shape(o)
    
    res%rx = r(0, s(1))
    res%ry = r(0, s(2))
    res%rz = r(0, s(3))    
  end function

  function array_lbox(o) result(res)
    implicit none
    type(array), intent(in) :: o
    type(box_info) :: res
    integer :: xs, ys, zs, xl, yl, zl

    call petsc_get_corners(o%data, res)
  end function

  function is_scalar_array(o) result(res)
    implicit none
    type(array), intent(in) :: o
    logical :: res
    integer :: s(3)

    s = shape(o)
    res = (s(1)==1 .and. &
         s(2)==1 .and. s(3)==1)
  end function

  function is_valid_array(t) result(res)
    implicit none
    type(array) :: t
    logical :: res
    
    if(any(shape(t) <= 0) &
         .or. t%data /=0) then
       res = .true.
    else
       res = .false.
    end if
  end function

  subroutine disp_info_array(o, prefix)
    implicit none
    type(array), intent(in) :: o
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
       write(*, "(4X, A, A)") "var_type = ", o%var_type
       write(*, "(4X, A, I3.1)") "ref count :", o%ref_cnt
#endif
       call disp_shape(o%m_shape, '4X')
       
       write(*, "(4X, A, L1)") "is_field    : ", o%is_field
       if(o%is_field) &
            write(*, "(4X, A, I0.2)") &
            "grid_pos : ", o%grid_pos
       write(*, "(4X, A, Z16.16)") "data : 0x", o%data
    end if
  end subroutine
  
end module 

