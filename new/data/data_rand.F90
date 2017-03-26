
subroutine data_rand(A, t_shape, ierr)
    implicit none
#include "petsc.h"

    Vec, intent(out):: A
    DM :: ada
    PetscErrorCode, intent(out) ::ierr
    integer :: t_shape(:)
    PetscInt :: m, n, k
    integer :: num_dim
    PetscInt :: xs, ys, zs, xl, yl, zl
    PetscInt::  ista,iend,ilocal
    PetscInt,   allocatable::idxm(:),idxn(:)
    PetscScalar,allocatable::row(:)
    integer :: i,j
    PetscLogEvent            ::  ievent
    integer :: dof = 1, s = 1
    PetscRandom :: rctx
    
    call PetscLogEventRegister("data_rand",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    num_dim = size(t_shape)

    m = t_shape(1)
    n = t_shape(2)
    k = t_shape(3)

    call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,         &
         DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,                      &
         DMDA_STENCIL_STAR, m,n,k,PETSC_DECIDE,PETSC_DECIDE,     &
         PETSC_DECIDE,dof,s,                                     &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                  &
         PETSC_NULL_INTEGER,ada,ierr)
   
    call DMCreateGlobalVector(ada, A, ierr)
    call PetscRandomCreate(PETSC_COMM_WORLD, rctx, ierr)    
    call VecSetRandom(A, rctx)
    call PetscLogEventEnd(ievent,ierr)
 end subroutine
