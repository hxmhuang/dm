program main 

    use dm 
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

    type(Matrix)    ::  A,B    
    PetscErrorCode  ::  ierr 
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    
    A=dm_ones(4,3)
    B = A

    ierr=dm_view(A)
    ierr=dm_view(B)

    ierr=dm_destroy(A)
    ierr=dm_destroy(B)

    call PetscFinalize(ierr)
end program
