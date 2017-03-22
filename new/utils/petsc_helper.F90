
#define I1 :
#define I2 :,:
#define I3 :,:,:

module petsc_helper
  use ot_geom

#:for d in [1,2,3]
  interface get_local_arr
     module procedure get_local_arr${d}$
  end interface get_local_arr
#:endfor

#:for d in [1,2,3]
  interface restore_local_arr
     module procedure restore_local_arr${d}$
  end interface restore_local_arr
#:endfor
  
contains
    subroutine vec_get_dim(A, num_dim)
    implicit none
#include "petsc.h"
    Vec :: A
    DM :: DM_A
    integer :: ierr
    integer, intent(out) :: num_dim
    
    call VecGetDM(A, DM_A, ierr)

    call DMDAGetInfo(DM_A, num_dim, &
         PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,ierr)
  end subroutine

  subroutine dm_get_dim(DM_A, num_dim)
    implicit none
#include "petsc.h"
    DM :: DM_A
    integer :: ierr
    integer, intent(out) :: num_dim
    
    call DMDAGetInfo(DM_A, num_dim, &
         PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,ierr)
  end subroutine

!   subroutine dm_get_corners(DM_A, xs, ys, zs, xl, yl, zl)
!     implicit none
! #include "petsc.h"
!     DM :: DM_A
!     integer :: ierr
!     integer, intent(out) :: xs, ys, zs, xl, yl, zl
    
!     call DMDAGetCorners(DM_A, xs, ys, zs, xl, yl, zl,ierr)

!   end subroutine

  subroutine dm_get_corners(A, b)
    implicit none
    DM :: A
    integer :: xs, ys, zs, xl, yl, zl
    type(box_info), intent(out) :: b
    integer :: ierr
    
    call DMDAGetCorners(A, xs, ys, zs, xl, yl, zl,ierr)    
    b%starts = (/xs, ys, zs/)
    b%ends = (/xs+xl-1, ys+yl-1, zs+zl-1/)
  end subroutine

  subroutine vec_get_corners(A, b)
    implicit none
    Vec :: A
    DM :: DM_A
    integer :: xs, ys, zs, xl, yl, zl
    type(box_info), intent(out) :: b
    integer :: ierr
    
    call VecGetDM(A, DM_A, ierr)
    call DMDAGetCorners(DM_A, xs, ys, zs, xl, yl, zl,ierr)

    b%starts = (/xs, ys, zs/)
    b%ends = (/xs+xl-1, ys+yl-1, zs+zl-1/)
  end subroutine

#:for d in [1,2,3]
  subroutine get_local_arr${d}$(data, arr)
    implicit none
#include "petsc.h"    
    Vec :: data
    DM :: data_dm
    PetscScalar, pointer, intent(out) :: arr(I${d}$)
    integer :: ierr
    
    call VecGetDM(data, data_dm, ierr)

    call DMDAVecGetArrayF90(data_dm, data, arr, ierr)
  end subroutine

  subroutine restore_local_arr${d}$(data, arr)
    implicit none
#include "petsc.h"
    Vec :: data
    DM  :: data_dm
    PetscScalar, pointer, intent(inout) :: arr(I${d}$)
    integer :: ierr
    
    call VecGetDM(data, data_dm, ierr)

    call DMDAVecRestoreArrayF90(data_dm, data, arr, ierr)
  end subroutine  
#:endfor
  
end module

#undef I1
#undef I2
#undef I3
