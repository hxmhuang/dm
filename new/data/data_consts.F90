  subroutine data_consts(A, alpha, t_shape, ierr)
    implicit none
#include "petsc.h"
    Vec, intent(out):: A
    DM :: ada
    PetscErrorCode, intent(out) ::ierr
    integer :: t_shape(:)
    PetscInt :: m, n, k
    integer :: num_dim
    PetscInt :: xs, ys, zs, xl, yl, zl
    PetscScalar,intent(in)::  alpha
    PetscInt::  ista,iend,ilocal
    PetscInt,allocatable::idxm(:),idxn(:)
    PetscScalar,allocatable::row(:)
    integer :: i,j
    PetscLogEvent            ::  ievent
    PetscScalar, pointer :: x3(:,:,:)

    call PetscLogEventRegister("data_consts",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    m = t_shape(1)
    n = t_shape(2)
    k = t_shape(3)
    
    call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,               &
         DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,                      &
         DMDA_STENCIL_STAR, m,n,k,PETSC_DECIDE,PETSC_DECIDE,     &
         PETSC_DECIDE, 1, 0,                                     &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                  &
         PETSC_NULL_INTEGER,ada,ierr)

    call DMGetGlobalVector(ada, A,ierr)
    call DMDAGetCorners(ada,xs,ys,zs, xl,yl,zl,ierr)
    call DMDAVecGetArrayF90(ada, A, x3,ierr)

    x3(xs:xs+xl-1, ys:ys+yl-1, zs:zs+zl-1) = alpha;

    call DMDAVecRestoreArrayF90(ada, A, x3, ierr)
       
    call PetscLogEventEnd(ievent,ierr)
  end subroutine data_consts
