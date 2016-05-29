program main 

    use dm 
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

    type(Matrix)    :: A,B,C,D,E,F,G,H 
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
    
    if(myrank==0) then 
       print *, "============Input paramenters============"
        print *, "m=",m,",n=",n,"ep=",ep,"debug=",debug
     endif 
	
    
 	if(myrank==0) print *, "==============Test dm_zeros==============="
    A=dm_zeros(m,n)
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
 	endif
    ierr=dm_destroy(A)


 	if(myrank==0) print *, "==============Test dm_ones================"
    A=dm_ones(m,n)
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
 	endif
    ierr=dm_destroy(A)

 	if(myrank==0) print *, "==============Test dm_seq================"
    A=dm_seqs(m,n)
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
 	endif
    ierr=dm_destroy(A)


    if(myrank==0) print *, "==============Test dm_eye==============="
    A=dm_eyes(m,m)	
    if(debug) then
        if(myrank==0) print *, ">A=(m,m)"
        ierr=dm_view(A)
 	endif
 	ierr=dm_destroy(A)

    A=dm_eyes(m,m*2)	
    if(debug) then
        if(myrank==0) print *, ">A=(m,2m)"
        ierr=dm_view(A)
 	endif
 	ierr=dm_destroy(A)

    A=dm_eyes(2*m,m)	
    if(debug) then
        if(myrank==0) print *, ">A=(2m,m)"
        ierr=dm_view(A)
 	endif
 	ierr=dm_destroy(A)


 	if(myrank==0) print *, "==============Test dm_copy==============="
 	A=dm_ones(m,m)
    B=A
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)


    if(myrank==0) print *, "==============Test dm_add==============="
 	A=dm_eyes(m,m)
    B=dm_ones(m,m)
    !C=dm_add1(A,B)
    C=A+B
    D=dm_eyes(m,m)+dm_ones(m,m)
    E=dm_eyes(m,m)+B
    F=A+dm_ones(m,m)
    G=A+A+A
    H=B+G
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C="
        ierr=dm_view(C)
        if(myrank==0) print *, ">D="
        ierr=dm_view(D)
        if(myrank==0) print *, ">E="
        ierr=dm_view(E)
        if(myrank==0) print *, ">F="
        ierr=dm_view(F)
        if(myrank==0) print *, ">G="
        ierr=dm_view(G)
        if(myrank==0) print *, ">H="
        ierr=dm_view(H)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
 	ierr=dm_destroy(F)
 	ierr=dm_destroy(G)
 	ierr=dm_destroy(H)








    call PetscFinalize(ierr)
end program
