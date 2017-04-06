  subroutine data_seqs(data, t_shape, ierr)
    implicit none
#include "petsc.h"
    Vec, intent(out):: data
    integer :: t_shape(3)
    DM :: data_dm
    PetscErrorCode, intent(out) ::ierr
    PetscInt :: m, n, k
    integer :: num_dim
    PetscInt :: xs, ys, zs, xl, yl, zl, xe, ye, ze, x, y, z
    PetscInt::  ista,iend,ilocal
    PetscInt,allocatable::idxm(:),idxn(:)
    PetscScalar,allocatable::row(:)
    integer :: i,j
    PetscLogEvent            ::  ievent
    integer :: dof = 1, s = 1
    PetscScalar, pointer :: x1(:), x2(:,:), x3(:,:,:)

    call PetscLogEventRegister("data_seqs",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    m = t_shape(1)
    n = t_shape(2)
    k = t_shape(3)

    call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,         &
         DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,                      &
         DMDA_STENCIL_STAR, m,n,k,PETSC_DECIDE,PETSC_DECIDE,     &
         PETSC_DECIDE,1, 0,                                      &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                  &
         PETSC_NULL_INTEGER,data_dm,ierr)

    call DMCreateGlobalVector(data_dm, data, ierr)
    call DMDAGetCorners(data_dm,xs,ys,zs,xl,yl,zl,ierr)
    call DMDAVecGetArrayF90(data_dm, data, x3, ierr)

    xe = xs + xl - 1
    ye = ys + yl - 1
    ze = zs + zl - 1

    ! print*, "m=", m, "n=", n, "k=",k
    ! print*, "xs = ", xs, "ys=", ys, "zs=", zs
    ! print*, "xe = ", xe, "ye=", ye, "ze=", ze
    
    !x3(zs:xs+xl-1, ys:ys+yl-1, zs:zs+zl-1) = alpha;

    do z = zs, ze
       do y = ys, ye
          do x = xs, xe
             x3(x, y, z) = x + y * m + z * m * n
          end do
       end do
    end do

    call DMDAVecRestoreArrayF90(data_dm, data, x3, ierr)
    call DMDestroy(data_dm, ierr)
    call PetscLogEventEnd(ievent,ierr)
  end subroutine
