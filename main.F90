program main 

    use dmc 
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
    
    A=mat_ones(4,3)
    B = A

    ierr=mat_view(A)
    ierr=mat_view(B)

    ierr=mat_destroy(A)
    ierr=mat_destroy(B)

    call PetscFinalize(ierr)
end program
