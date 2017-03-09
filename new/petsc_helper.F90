module petsc_helper
  
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

  subroutine dm_get_corners(DM_A, xs, ys, zs, xl, yl, zl)
    implicit none
#include "petsc.h"
    DM :: DM_A
    integer :: ierr
    integer, intent(out) :: xs, ys, zs, xl, yl, zl
    
    call DMDAGetCorners(DM_A, xs, ys, zs, xl, yl, zl,ierr)

  end subroutine

  subroutine vec_get_corners(A, xs, ys, zs, xl, yl, zl)
    implicit none
#include "petsc.h"
    Vec :: A
    DM :: DM_A
    integer :: ierr
    integer, intent(out) :: xs, ys, zs, xl, yl, zl

    call VecGetDM(A, DM_A, ierr)
    
    call DMDAGetCorners(DM_A, xs, ys, zs, xl, yl, zl,ierr)

  end subroutine
  
end module
