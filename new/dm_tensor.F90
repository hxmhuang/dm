
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 

module dm_tensor
  use dm_common
  use dm_data

  public :: tensor_new, tensor_destroy
  
  integer, parameter :: zero3(3) = (/0,0,0/)
  
  type tensor
     !the left and right leaf for the tree
     type(tensor), pointer :: left => null(), right => null()

     !node type, 0 for data, 1 for '+', 2 for '-', 3 for '*', 4 for '/'
     integer :: node_type = 0

     !dim and shape should not be specified by user!
     integer :: m_dim = 0, m_shape(3) = zero3 

     !the pointers used to point petsc raw data
     real(kind=8), pointer :: data1d(:) => null(), &
          data2d(:,:) => null(), data3d(:,:,:) => null()
     
     Vec :: data = 0
     !DM  :: data_dm = 0

     !check if it is an implicit variable
     logical :: is_implicit = .false.
  end type tensor
  
contains

  subroutine display(objA, prefix)
    type(tensor), intent(in),target :: objA
    type(tensor), pointer :: A
    character(len=*), optional, intent(in) :: prefix
    real(kind=8), pointer :: x1(:), x2(:,:), x3(:,:,:)
    integer :: i

    A => objA

    write(*,*) ""
    if(present(prefix)) then
       write(*, "(A)") prefix
    else
       write(*, "(A)") "ans = "
    endif
    
    write(*, "(4X, A, 100I2)") "node_type : ", A%node_type
    
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
       write(*, "(4X, A, Z16.6)"), "left  = 0x", loc(A%left)
       write(*, "(4X, A, Z16.6)"), "right = 0x", loc(A%right)       
    end if
    write(*,*) ""    
  end subroutine display

  function tensor_new(m_shape) result(A)
    implicit none
    integer, intent(in) :: m_shape(:)
    type(tensor):: A

    A%m_dim = size(m_shape)
    A%m_shape(1:A%m_dim) = m_shape
  end function

  subroutine tensor_destroy(A, ierr)
    implicit none
    type(tensor), intent(inout) :: A
    integer, intent(out) :: ierr

    call data_destroy(A%data, ierr)
  end subroutine 
  
  subroutine tensor_plus(A, B, C) !A = B + C
    implicit none
    type(tensor), intent(inout) :: A
    type(tensor), intent(inout) :: B, C
    integer :: ierr

    if(B%is_implicit) then
       call tensor_copy(A, B)
       call data_plus(A%data, C%data, ierr) !A = A + C
       if(C%is_implicit) then
          call tensor_destroy(C, ierr)
       endif
    else
       if(C%is_implicit) then
          call tensor_copy(A, C)
          call data_plus(A%data, B%data, ierr)
       else
          call tensor_duplicate(A, B)
          call data_plus(A%data, B%data, C%data, ierr)
       end if
    end if
    
  end subroutine 

  !A = B - C
  subroutine tensor_minus(A, B, C)
    implicit none
    type(tensor), intent(inout) :: A
    type(tensor), intent(inout) :: B, C
    integer :: ierr

    if(B%is_implicit) then
       call tensor_copy(A, B)
       call data_minus(A%data, C%data, ierr) !A = A + C
       if(C%is_implicit) then
          call tensor_destroy(C, ierr)
       endif
    else
       if(C%is_implicit) then
          call tensor_copy(A, C)
          call data_minus(A%data, B%data, ierr)
       else
          call tensor_duplicate(A, B)
          call data_minus(A%data, B%data, C%data, ierr)
       end if
    end if
    
  end subroutine 

  subroutine tensor_mult(A, B, C)
    implicit none
    type(tensor), intent(inout) :: C
    type(tensor), intent(inout) :: A, B
    integer :: ierr
    
    if(B%is_implicit) then
       call tensor_copy(A, B)
       call data_mult(A%data, C%data, ierr) !A = A + C
       if(C%is_implicit) then
          call tensor_destroy(C, ierr)
       endif
    else
       if(C%is_implicit) then
          call tensor_copy(A, C)
          call data_mult(A%data, B%data, ierr)
       else
          call tensor_duplicate(A, B)
          call data_mult(A%data, B%data, C%data, ierr)
       end if
    end if
  end subroutine 

  subroutine tensor_divd(A, B, C)
    implicit none
    type(tensor), intent(inout) :: C
    type(tensor), intent(inout) :: A, B
    integer :: ierr
    
    if(B%is_implicit) then
       call tensor_copy(A, B)
       call data_divd(A%data, C%data, ierr) !A = A + C
       if(C%is_implicit) then
          call tensor_destroy(C, ierr)
       endif
    else
       if(C%is_implicit) then
          call tensor_copy(A, C)
          call data_divd(A%data, B%data, ierr)
       else
          call tensor_duplicate(A, B)
          call data_divd(A%data, B%data, C%data, ierr)
       end if
    end if
  end subroutine 

  !naively copy the data member and pointers from B to A
  subroutine tensor_copy(A, B)
    implicit none
    type(tensor), intent(inout) :: A
    type(tensor), intent(in) :: B
    integer :: ierr
    
    A%m_dim = B%m_dim
    A%m_shape = B%m_shape
    A%left  => B%left
    A%right => B%right
    A%node_type = B%node_type
    
    call data_destroy(A%data,  ierr)

    A%data = B%data
  end subroutine

  !deep copy tensor, all data are copied from B to A
  subroutine tensor_deep_copy(A, B)
    implicit none
    type(tensor), intent(out) :: A
    type(tensor), intent(in) :: B
    integer :: ierr
    
    A%m_dim = B%m_dim
    A%m_shape = B%m_shape
    A%left => B%left
    A%right => B%right
    A%node_type = B%node_type
    
    call data_destroy(A%data,  ierr) !safe destroy
    call data_clone(A%data, B%data, ierr)
  end subroutine tensor_deep_copy

  !duplicate tensor structure, but the data is not copied
  subroutine tensor_duplicate(A, B)
    implicit none
    type(tensor), intent(out) :: A
    type(tensor), intent(in) :: B
    integer :: ierr
    
    A%m_dim = B%m_dim
    A%m_shape = B%m_shape
    A%left => B%left
    A%right => B%right
    A%node_type = B%node_type
    
    call data_destroy(A%data,  ierr) !safe destroy
    call data_duplicate(A%data,  B%data,  ierr)
  end subroutine 
  
end module 

