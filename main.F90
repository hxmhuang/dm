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
    debug=.false.

    ierr=dm_init()
    
    myrank=dm_comm_rank()
    
    mysize=dm_comm_size()
    
    m=dm_get_int('-m')
    n=dm_get_int('-n')
    ep=dm_get_real('-ep')
    !debug=dm_get_bool('-debug')

    call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-debug',debug,PETSC_NULL_BOOL,ierr)
    
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
        if(myrank==0) print *, ">C=A+B"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_eyes(m,m)+dm_ones(m,m)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_eyes(m,m)+B"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=A+dm_ones(m,m)"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=A+A+A"
        ierr=dm_view(G)
        if(myrank==0) print *, ">H=B+G"
        ierr=dm_view(H)
 	endif
    !TODO:There is a bug to free matrix A.
    !A=A+A	
    !if(debug) then
    !    if(myrank==0) print *, ">A="
    !    ierr=dm_view(a)
 	!endif
    ierr=dm_destroy(A)
  	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
 	ierr=dm_destroy(F)
 	ierr=dm_destroy(G)
 	ierr=dm_destroy(H)


    if(myrank==0) print *, "==============Test dm_hjoin==============="
  	A=dm_eyes(m,m)
    B=dm_ones(m,m)
    C=A .hjoin. B
    D=dm_eyes(m,m) .hjoin. dm_ones(m,m)
    E=dm_eyes(m,m) .hjoin. B
    F=A .hjoin. dm_ones(m,m)
    G=A .hjoin. A .hjoin. A
    H=B .hjoin. G
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=A .hjoin. B"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_eyes(m,m) .hjoin. dm_ones(m,m)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_eyes(m,m) .hjoin. B"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=A .hjoin. dm_ones(m,m)"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=A .hjoin. A .hjoin. A"
        ierr=dm_view(G)
        if(myrank==0) print *, ">H=B .hjoin. G"
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


    if(myrank==0) print *, "==============Test dm_mult==============="
    A=dm_ones(m,m)
    B=dm_eyes(m,m*2)
    C=A*B
    D=A*(dm_eyes(m,m*2))
    E=dm_ones(m,m)*B
    F=dm_ones(m,m)*dm_eyes(m,m*2) 
    G=A*A
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=A*B"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=A*dm_eyes(m,m*2)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_ones(m,m)*B"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=dm_ones(m,m)*dm_eyes(m,m*2)"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=A*A"
        ierr=dm_view(G)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
    ierr=dm_destroy(F)
 	ierr=dm_destroy(G)


    if(myrank==0) print *, "==============Test dm_eprod==============="
    A=dm_seqs(m,m)
    B=dm_eyes(m,m)
    C=A .eprod. B
    D=A .eprod. (dm_eyes(m,m))
    E=dm_seqs(m,m) .eprod. B
    F=dm_seqs(m,m) .eprod. dm_eyes(m,m) 
    G=A .eprod. A
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=A.*B"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=A.*dm_eyes(m,m*2)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_seqs(m,m).*B"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=dm_seqs(m,m).*dm_eyes(m,m*2)"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=A.*A"
        ierr=dm_view(G)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
    ierr=dm_destroy(F)
 	ierr=dm_destroy(G)


    if(myrank==0) print *, "==============Test dm_rep==============="
    A=dm_eyes(m,m)
    B=dm_rep(A,3,2) 
    C=dm_rep(dm_eyes(m,m),3,2) 
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B=dm_rep(A,3,2)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_rep(dm_eyes(m,m),3,2)"
        ierr=dm_view(C)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_sum==============="
    A=dm_seqs(m,n)
    B=dm_sum(A,1)
    C=dm_sum(A,2)
    D=dm_sum(dm_seqs(m,n),1)
    E=dm_sum(dm_seqs(m,n),2)
    if(debug) then
        if(myrank==0) print *, ">A=dm_seqs(m,n)"
        ierr=dm_view(A)
        if(myrank==0) print *, ">B=dm_sum(A,1)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_sum(A,2)"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_sum(dm_seqs(m,n),2)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_sum(dm_seqs(m,n),2)"
        ierr=dm_view(E)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)



    call PetscFinalize(ierr)
end program
