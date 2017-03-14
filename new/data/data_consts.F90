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
    integer :: dof = 1, s = 1

    PetscScalar, pointer :: x1(:), x2(:,:), x3(:,:,:)

    call PetscLogEventRegister("data_consts",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    num_dim = size(t_shape)

    m = t_shape(1)
    if(num_dim >= 2) n = t_shape(2)
    if(num_dim >= 3) k = t_shape(3)

    if(num_dim == 1) then
       call DMDACreate1d(PETSC_COMM_WORLD, &
         &     DM_BOUNDARY_NONE,        &
         &     m,  dof, s,     &
         &     PETSC_NULL_INTEGER,ada,ierr)

       
       call DMGetGlobalVector(ada, A,ierr)
       call DMDAGetCorners(ada,xs, PETSC_NULL_INTEGER, &
            PETSC_NULL_INTEGER, xl, &
            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)       
       call DMDAVecGetArrayF90(ada, A, x1,ierr)

       !print*, "xs=", xs, "xl=", xl
       x1(xs:xs+xl-1) = alpha;
       
       call DMDAVecRestoreArrayF90(ada, A, x1, ierr)
       !call VecView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
       !call DMRestoreGlobalVector(ada, A, ierr)
       !call DMDestroy(ada,ierr)
       
    else if (num_dim == 2) then

       call DMDACreate2d(PETSC_COMM_WORLD, &
         &     DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,        &
         &     DMDA_STENCIL_STAR, m,n, PETSC_DECIDE,      &
         &     PETSC_DECIDE,dof,s,                       &
         &     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,    &
         &     ada,ierr)
       
       call DMGetGlobalVector(ada, A,ierr)
       call DMDAGetCorners(ada,xs, ys, &
            PETSC_NULL_INTEGER, xl, yl, &
            PETSC_NULL_INTEGER, ierr)       
       call DMDAVecGetArrayF90(ada, A, x2,ierr)
       
       
       x2(xs:xs+xl-1, ys:ys+yl-1) = alpha

       call DMDAVecRestoreArrayF90(ada, A, x2, ierr)
       !call VecView(A, PETSC_VIEWER_STDOUT_WORLD,ierr)
       !call DMRestoreGlobalVector(ada, A, ierr)
       !call DMDestroy(ada,ierr)

    else if(num_dim == 3) then
       
       call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,               &
            &     DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,                      &
            &     DMDA_STENCIL_STAR, m,n,k,PETSC_DECIDE,PETSC_DECIDE,     &
            &     PETSC_DECIDE,dof,s,                                     &
            &     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                  &
            &     PETSC_NULL_INTEGER,ada,ierr)

       call DMGetGlobalVector(ada, A,ierr)
       call DMDAGetCorners(ada,xs,ys,zs, xl,yl,zl,ierr)
       call DMDAVecGetArrayF90(ada, A, x3,ierr)

       x3(zs:xs+xl-1, ys:ys+yl-1, zs:zs+zl-1) = alpha;
       
       call DMDAVecRestoreArrayF90(ada, A, x3, ierr)
       !call VecView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
       !call DMRestoreGlobalVector(ada, A, ierr)
       !call DMDestroy(ada,ierr)
    endif

    call PetscLogEventEnd(ievent,ierr)
  end subroutine data_consts
