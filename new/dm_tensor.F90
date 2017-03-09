
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 
#:include "type_def.fypp"

#include "dm_type.h"

module dm_tensor
  use dm_common
  use dm_data
  use func_helper
  
  public :: tensor_new, tensor_destroy
  
  integer, parameter :: zero3(3) = (/0,0,0/)

  type grid
     character(1) :: grid_type !'B' or 'C'
     type(tensor), pointer :: dx_3d => null()
     type(tensor), pointer :: dy_3d => null()
     type(tensor), pointer :: dz_3d => null()     
  end type grid
  
  type tensor
     !dim and shape should not be specified by user!
     integer :: m_dim = 0, m_shape(3) = zero3 
     
     Vec :: data = 0

     !check if it is an implicit variable
     logical :: is_implicit = .false.

     logical :: is_field = .false.
     integer :: grid_pos

  end type tensor

  type tensor_ptr
     type(tensor), pointer :: ptr
  end type tensor_ptr
  
  !here we use a trick to extract range information from an empty array
  !see find_range function
  integer, allocatable :: r(:)

  interface
     subroutine reg_func(p, id) &
          bind(C, name="reg_func")
       implicit none
       C_POINTER :: p
       integer :: id
     end subroutine
  end interface
  
contains

  subroutine init_tensor(ierr)
    implicit none
    integer, intent(out) :: ierr

    !initialize data module
    call init_data(ierr)

#:for e in L
#:if e[1] >= 50
    !write(*, "(A, Z16.6)"), "f_addr=", loc(tensor_${e[2]}$)
    call reg_func(loc(tensor_${e[2]}$), ${e[0]}$)
#:endif
#:endfor

  end subroutine
  
  function find_range(range) result(res)
    implicit none
    integer, intent(inout) :: range(:)
    integer :: res(2)

    !print*, size(range)
    res(1) = (loc(range) - loc(r)) / 4
    res(2) = ((loc(range) + size(range) * 4) - loc(r) - 4)/4

  end function
  
  subroutine display2(objA, prefix)
#include "petsc.h"            
    type(tensor), intent(in),target :: objA
    type(tensor), pointer :: A
    character(len=*), optional, intent(in) :: prefix
    real(kind=8), pointer :: x1(:), x2(:,:), x3(:,:,:)
    integer :: i

    A => objA

    print*, "adsafasfasf"
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
  subroutine tensor_${op[2]}$(res, operands, alpha, beta, args) 
    implicit none
    
    type(tensor), intent(inout), pointer :: res
    type(tensor_ptr), intent(in)  :: operands(:)
    Vec,allocatable :: operands_vec(:)
    real(8), intent(in) :: alpha(:), beta(:)
    real(8), intent(in) :: args(10)
    integer :: n
    integer :: ierr, i
    
    allocate(operands_vec(size(operands)))

    if(.not. associated(res)) then
       allocate(res)
    endif
    
    call tensor_duplicate(res, operands(1)%ptr)
    
    do i = 1, size(operands)
       operands_vec(i) = operands(i)%ptr%data
    enddo
    
    call data_${op[2]}$(res%data, operands_vec, alpha, beta, args, ierr)


    ! if(B%is_implicit) then
    !    call tensor_copy(A, B)
    !    call data_plus(A%data, C%data, ierr) !A = A + C
    !    if(C%is_implicit) then
    !       call tensor_destroy(C, ierr)
    !    endif
    ! else
    !    if(C%is_implicit) then
    !       call tensor_copy(A, C)
    !       call data_plus(A%data, B%data, ierr)
    !    else
    !       call tensor_duplicate(A, B)
    !       call data_plus(A%data, B%data, C%data, ierr)
    !    end if
    ! end if
  end subroutine
  #:endif
#:endfor

  subroutine tensor_copy_structure(A, B)
    implicit none
    type(tensor), intent(inout) :: A
    type(tensor), intent(in) :: B

    A%m_dim   = B%m_dim
    A%m_shape = B%m_shape
    
  end subroutine
  
  !> naively copy the data member and pointers from B to A
  subroutine tensor_copy(A, B)
    implicit none
    type(tensor), intent(inout) :: A
    type(tensor), intent(in) :: B
    integer :: ierr
    
    call tensor_copy_structure(A, B)
    call data_destroy(A%data,  ierr)

    A%data = B%data
  end subroutine

  !> deep copy tensor, all data are copied from B to A
  subroutine tensor_deep_copy(A, B)
    implicit none
    type(tensor), intent(out) :: A
    type(tensor), intent(in)  :: B
    integer :: ierr

    call tensor_copy_structure(A, B)
    call data_destroy(A%data,  ierr) !safe destroy
    call data_clone(A%data, B%data, ierr)
  end subroutine tensor_deep_copy

  !> duplicate tensor structure, the data is allocated but not copied
  subroutine tensor_duplicate(A, B)
    implicit none
    type(tensor), intent(out) :: A
    type(tensor), intent(in) :: B
    integer :: ierr
    
    call tensor_copy_structure(A, B)
    call data_destroy(A%data,  ierr) !safe destroy
    call data_duplicate(A%data,  B%data,  ierr)
  end subroutine 

  
end module 

