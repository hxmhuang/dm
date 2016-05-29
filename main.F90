program main 

    use dm 
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

    type(Matrix)    :: A,B    
    integer         :: myrank, mysize 
    integer         :: m,n 
    real(kind=8)    :: ep
    logical         :: debug 
    integer         :: ierr

    ierr=dm_init()
    
    myrank=dm_comm_rank()
    
    mysize=dm_comm_size()
    
    m=dm_get_int('-m')
    n=dm_get_int('-n')
    ep=dm_get_real('-ep')
    debug=dm_get_bool('-debug')
    print *,m,n,ep,debug
    
    A=dm_ones(4,3)
    
    B = A

    ierr=dm_view(A)
    ierr=dm_view(B)

    ierr=dm_destroy(A)
    ierr=dm_destroy(B)

    call PetscFinalize(ierr)
end program
