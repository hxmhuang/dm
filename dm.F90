! -----------------------------------------------------------------------
! Distributed Matrix Computing Library
! -----------------------------------------------------------------------
module dm 
    use dm_type
    implicit none
    interface assignment(=)
        module procedure dm_copyEx
        module procedure dm_copyIm 
    end interface

    interface dm_destroy
        module procedure dm_destroyEx 
        module procedure dm_destroyIm 
    end interface

contains

! -----------------------------------------------------------------------
!Create a matrix with m*n size
! -----------------------------------------------------------------------
function dm_create(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n 
    type(MatrixIm)                ::  A
    PetscErrorCode              ::  ierr
    ! generate matrix A with size m*n
    call MatCreate(PETSC_COMM_WORLD,A%x,ierr);
    call MatSetSizes(A%x,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr)
    call MatSetFromOptions(A%x,ierr)
    call MatSetUp(A%x,ierr)
    call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)
end function 

! -----------------------------------------------------------------------
!Destroy a matrix to free the memory
! -----------------------------------------------------------------------
function dm_destroyEx(A) result(ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),   intent(in)  ::  A
    PetscErrorCode              ::  ierr
    ! destroy matrix A
    call MatDestroy(A%x,ierr)
end function 

! -----------------------------------------------------------------------
!Destroy a implict matrix to free the memory
! -----------------------------------------------------------------------
function dm_destroyIm(A) result(ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(MatrixIm),   intent(in)  ::  A
    PetscErrorCode              ::  ierr
    ! destroy matrix A
    call MatDestroy(A%x,ierr)
end function 

! -----------------------------------------------------------------------
! Print a matrix on screen
! -----------------------------------------------------------------------
function dm_view(A) result(ierr) 
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),  intent(in)   ::  A
    PetscErrorCode              ::  ierr
    ! veiw matrix A
    call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)

    call MatView(A%x,PETSC_VIEWER_STDOUT_WORLD, ierr)
end function 

! -----------------------------------------------------------------------
! A=0 
! -----------------------------------------------------------------------
function dm_zeros(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n 
    type(MatrixIm)                ::  A
    PetscErrorCode              ::  ierr

    !call mat_create(A,m,n,ierr)
    A=dm_create(m,n)
    call MatZeroEntries(A%x,ierr)
end function


function dm_ones(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n
    PetscErrorCode              ::  ierr
    type(MatrixIm)                ::  A

    PetscInt                    ::  ista,iend
    PetscInt,allocatable        ::  idxn(:)
    PetscScalar,allocatable     ::  row(:),results(:)
    integer                     ::  i,j

    A=dm_create(m,n)
    call MatGetOwnershipRange(A%x,ista,iend,ierr)
    !print *,">ista=",ista,"iend=",iend,">ncol=",ncol
    allocate(idxn(n),row(n),results(n))

    !call MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY,ierr)
    !call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)
    do i=ista,iend-1
        do j=1,n
            idxn(j)=j-1
            row(j)=1.0
        enddo
        call MatSetValues(A%x,1,i,n,idxn,row,INSERT_VALUES,ierr)
    enddo
    call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)

    deallocate(idxn,row,results)
end function

! -----------------------------------------------------------------------
! B=A. This function uses the implicit matrix A directly because A is not need to free. 
! -----------------------------------------------------------------------
subroutine dm_copyIm(B,A)
    implicit none

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),            intent(out) ::  B
    type(MatrixIm),            intent(in)  ::  A
    PetscErrorCode              ::  ierr
    B%x=A%x
end subroutine


subroutine dm_copyEx(B,A)
    implicit none

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),            intent(out) ::  B
    type(Matrix),            intent(in)  ::  A
    PetscErrorCode              ::  ierr
    
    call MatDuplicate(A%x,MAT_COPY_VALUES,B%x,ierr)
end subroutine


end module 


