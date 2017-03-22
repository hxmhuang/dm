#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 
#:include "type_def.fypp"
#include "node_type.h"

module ot_tensor
  use ot_common
  use ot_data
  use ot_type  
  use ot_ref

  public :: tensor_new, tensor_destroy
  
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

  interface dim
     module procedure tensor_dim
  end interface dim

  interface box
     module procedure tensor_box
  end interface box

  interface lbox
     module procedure tensor_lbox
  end interface lbox

  interface disp_info
     module procedure disp_info_tensor
  end interface disp_info

contains

  subroutine init_tensor(ierr)
    implicit none
    integer, intent(out) :: ierr

    !initialize data module
    call init_data(ierr)

  end subroutine
  
  function find_range(range) result(res)
    implicit none
    integer, intent(inout) :: range(:)
    integer :: res(2)

    !print*, size(range)
    res(1) = (loc(range) - loc(r)) / 4
    res(2) = ((loc(range) + size(range) * 4) - loc(r) - 4)/4

  end function

  function tensor_dim(o) result(res)
    implicit none
    type(tensor), intent(in) :: o
    integer :: res
    res = o%m_dim
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

    res%starts = 0
    res%ends = shape(o) - 1
  end function

  function tensor_lbox(o) result(res)
    implicit none
    type(tensor), intent(in) :: o
    type(box_info) :: res
    integer :: xs, ys, zs, xl, yl, zl

    call vec_get_corners(o%data, res)

  end function

  ! subroutine tensor_ensure_valid(o)
  !   implicit none
  !   type(tensor), intent(in) :: o

  !   call assert(o%data /= 0, __FILE__, __LINE__, &
  !        "tensor must be associated with data")
    
  ! end subroutine
  
  subroutine disp_info_tensor(o, prefix)
    implicit none
    type(tensor), intent(in) :: o
    character(len=*), optional, intent(in) :: prefix
    integer :: i
    
    if(get_rank() == 0) then
       write(*,*) ""
       if(present(prefix)) then
          write(*, "(A)") prefix
       else
          write(*, "(A)") "ans = "
       endif

       write(*, "(4X, A, Z16.16)") "obj addr : 0x", loc(o)

       write(*, "(4X, A)", advance="no") "shape : ["
       do i = 1, o%m_dim
          write(*, "(I0.1)", advance="no") o%m_shape(i)
          if(i < o%m_dim) write(*, "(A)", advance="no") "x"
       enddo
       write(*, "(A)") "]"

       write(*, "(4X, A, L1)") "is_implicit : ", o%is_implicit
       write(*, "(4X, A, L1)") "is_field    : ", o%is_field
       if(o%is_field) write(*, "(4X, A, I0.2)") "grid_pos : ", o%grid_pos
       write(*, "(4X, A, Z16.16)") "data : 0x", o%data
    end if
  end subroutine
  
  subroutine display2(objA, prefix)
#include "petsc.h"            
    type(tensor), intent(in),target :: objA
    type(tensor), pointer :: A
    character(len=*), optional, intent(in) :: prefix
    real(kind=8), pointer :: x1(:), x2(:,:), x3(:,:,:)
    integer :: i

    A => objA

    !call VecView(objA%data, PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    write(*,*) ""
    if(present(prefix)) then
       write(*, "(A)") prefix
    else
       write(*, "(A)") "ans = "
    endif
    
    !write(*, "(4X, A, 100I2)") "node_type : ", A%node_type
    
    write(*, "(4X, A)", advance="no") "shape : ["
    do i = 1, A%m_dim
       write(*, "(I0.1)", advance="no") A%m_shape(i)
       if(i < A%m_dim) write(*, "(A)", advance="no") "x"
    enddo
    write(*, "(A)") "]"
    
    write(*, "(4X, A, L1)") "is_implicit : ", A%is_implicit
    write(*, "(A)") ""

    if(A%data /= 0) then
       select case (A%m_dim)
       case (1)
          call get_local_ptr(A%data,  x1)
          call dm_print("data=", x1)
       case (2)
          call get_local_ptr(A%data,  x2)
          call dm_print("data=", x2)          
       case (3)
          call get_local_ptr(A%data,  x3)
          call dm_print("data=", x3)
       end select
    else
       !write(*, "(4X, A, Z16.6)"), "left  = 0x", loc(A%left)
       !write(*, "(4X, A, Z16.6)"), "right = 0x", loc(A%right)       
    end if
    write(*,*) ""    
  end subroutine

  subroutine tensor_new(A, m_shape)
    implicit none
    integer, intent(in) :: m_shape(:)
    type(tensor), intent(out) :: A

    A%m_dim = size(m_shape)
    A%m_shape(1:A%m_dim) = m_shape
  end subroutine

  subroutine tensor_destroy(A, ierr)
    implicit none
    type(tensor), intent(inout) :: A
    integer, intent(out) :: ierr

    call data_destroy(A%data, ierr)
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

    dst%m_dim   = src%m_dim
    dst%m_shape = src%m_shape
    dst%is_implicit = src%is_implicit
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

  subroutine bind_data(o, data)
    implicit none
    type(tensor), intent(inout) :: o
    Vec, intent(in) :: data
    integer :: d_shape(3)
    integer :: xs, xl, ys, yl, zs, zl

    o%data = data    
    call vec_get_corners(data, o%local_block)

  end subroutine
  
end module 

