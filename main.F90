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
    type(Matrix)    :: X,Y,Z,U 
    integer         :: myrank, mysize 
    integer         :: m,n
	integer			:: idxm(2),idxn(2)
    real(kind=8)    :: ep,alpha
    real(kind=8)    :: a1,a2,a3
    logical         :: debug 
    integer         :: ierr
    character(len=50):: filename
	real(kind=8),allocatable :: array(:)
    debug=.false.

    call dm_init(ierr)
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
	call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
	
	if(myrank==0) then 
       print *, "==============Input paramenters==========="
        print *, "m=",m,",n=",n,",ep=",ep,",debug=",debug
     endif 
	
    
 	if(myrank==0) print *, "==============Test dm_zeros==============="
    A=dm_zeros(m,n)
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
 	endif
    call dm_destroy(A,ierr)


 	if(myrank==0) print *, "==============Test dm_ones================"
    A=dm_ones(m,n)
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
 	endif
    call dm_destroy(A,ierr)


 	if(myrank==0) print *, "==============Test dm_constants==========="
    A=dm_constants(m,n,3)
    B=dm_constants(m,n,4.5)
    C=dm_constants(m,n,real(5.5,kind=8))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_constants(m,n,4.5)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_constants(m,n,real(4.5,kind=8))"
        call dm_view(C,ierr)
 	endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)


 	if(myrank==0) print *, "==============Test dm_seqs================"
    A=dm_seqs(m,n)
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
 	endif
    call dm_destroy(A,ierr)


 	if(myrank==0) print *, "==============Test dm_m2n================="
    A=dm_m2n(0,2)
    if(debug) then
        if(myrank==0) print *, ">A=dm_m2n(0,2)"
        call dm_view(A,ierr)
 	endif
    call dm_destroy(A,ierr)


    if(myrank==0) print *, "==============Test dm_eye================="
    A=dm_eyes(m,m)	
    if(debug) then
        if(myrank==0) print *, ">A=(m,m)"
        call dm_view(A,ierr)
 	endif
 	call dm_destroy(A,ierr)

    A=dm_eyes(m,m*2)	
    if(debug) then
        if(myrank==0) print *, ">A=(m,2m)"
        call dm_view(A,ierr)
 	endif
 	call dm_destroy(A,ierr)

    A=dm_eyes(2*m,m)	
    if(debug) then
        if(myrank==0) print *, ">A=(2m,m)"
        call dm_view(A,ierr)
 	endif
 	call dm_destroy(A,ierr)


 	if(myrank==0) print *, "==============Test dm_copy================"
 	A=dm_eyes(m,m)
    B=A
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)


    if(myrank==0) print *, "==============Test dm_add================="
  	A=dm_eyes(m,m)
    B=dm_eyes(m,m)
    C=A+B
    D=dm_eyes(m,m)+dm_eyes(m,m)
    E=dm_eyes(m,m)+B
    F=A+dm_eyes(m,m)
    G=A+A+A
    H=B+G
    X=A+2.0
    Y=2+A
    Z=real(2,8)+A
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A+B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_eyes(m,m)+dm_eyes(m,m)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eyes(m,m)+B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=A+dm_eyes(m,m)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A+A+A"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=B+G"
        call dm_view(H,ierr)
        if(myrank==0) print *, ">X=A+2.0"
        call dm_view(X,ierr)
        if(myrank==0) print *, ">Y=2+A"
        call dm_view(Y,ierr)
        if(myrank==0) print *, ">Z=real(2,8)+A"
        call dm_view(Z,ierr)
 	endif
     A=A+A	
     if(debug) then
         if(myrank==0) print *, ">A=A+A"
        call dm_view(A,ierr)
 	 endif
    call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
 	call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)
 	call dm_destroy(H,ierr)
 	call dm_destroy(X,ierr)
 	call dm_destroy(Y,ierr)
 	call dm_destroy(Z,ierr)


    if(myrank==0) print *, "==============Test dm_minus==============="
  	A=dm_ones(m,m)
    B=dm_eyes(m,m)
    C=A-B
    D=dm_eyes(m,m)-dm_eyes(m,m)
    E=dm_eyes(m,m)-B
    F=A-dm_eyes(m,m)
    G=A-A-A
    H=B-G
    X=9*A-2.0
    Y=2-A
    Z=real(2,8)-A
    U=(0-A)+2
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A-B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_eyes(m,m)-dm_eyes(m,m)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eyes(m,m)-B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=A-dm_eyes(m,m)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A-A-A"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=B-G"
        call dm_view(H,ierr)
        if(myrank==0) print *, ">X=9*A-2.0"
        call dm_view(X,ierr)
        if(myrank==0) print *, ">Y=2-A"
        call dm_view(Y,ierr)
        if(myrank==0) print *, ">Z=real(2,8)-A"
        call dm_view(Z,ierr)
        if(myrank==0) print *, ">U=(0-A)+2"
        call dm_view(U,ierr)
 	endif
    A=A-A	
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
 	endif
    call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
 	call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)
 	call dm_destroy(H,ierr)
 	call dm_destroy(X,ierr)
 	call dm_destroy(Y,ierr)
 	call dm_destroy(Z,ierr)
 	call dm_destroy(U,ierr)


    if(myrank==0) print *, "==============Test dm_hjoin==============="
  	A=dm_eyes(m,m)
    B=dm_eyes(m,m)
    C=A .hj. B
    D=dm_eyes(m,m) .hj. dm_eyes(m,m)
    E=dm_eyes(m,m) .hj. B
    F=A .hj. dm_eyes(m,m)
    G=A .hj. A .hj. A
    H=B .hj. G
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A .hjoin. B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_eyes(m,m) .hjoin. dm_eyes(m,m)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eyes(m,m) .hjoin. B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=A .hjoin. dm_eyes(m,m)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A .hjoin. A .hjoin. A"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=B .hjoin. G"
        call dm_view(H,ierr)
 	endif
    call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
 	call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)
 	call dm_destroy(H,ierr)


    if(myrank==0) print *, "==============Test dm_vjoin==============="
  	A=dm_eyes(m,m)
    B=dm_eyes(m,m)
    C=A .vj. B
    D=dm_eyes(m,m) .vj. dm_eyes(m,m)
    E=dm_eyes(m,m) .vj. B
    F=A .vj. dm_eyes(m,m)
    G=A .vj. A .vj. A
    H=B .vj. G
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A .vj. B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_eyes(m,m) .vj. dm_eyes(m,m)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eyes(m,m) .vj. B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=A .vj. dm_eyes(m,m)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A .vj. A .vj. A"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=B .vj. G"
        call dm_view(H,ierr)
 	endif
    call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
 	call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)
 	call dm_destroy(H,ierr)



    if(myrank==0) print *, "==============Test dm_mult================"
    A=dm_eyes(m,m)
    B=dm_eyes(m,m*2)
    C=A*B
    D=A*(dm_eyes(m,m*2))
    E=dm_eyes(m,m)*B
    F=dm_eyes(m,m)*dm_eyes(m,m*2) 
    G=A*A
    X=A*2.0
    Y=2*A
    alpha=3.0
    Z=alpha*A
    U=A*alpha
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A*B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=A*dm_eyes(m,m*2)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eyes(m,m)*B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=dm_eyes(m,m)*dm_eyes(m,m*2)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A*A"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">X=A*2.0"
        call dm_view(X,ierr)
        if(myrank==0) print *, ">Y=2*A"
        call dm_view(Y,ierr)
        if(myrank==0) print *, ">Z=alpha*A"
        call dm_view(Y,ierr)
        if(myrank==0) print *, ">U=A*alpha"
        call dm_view(Y,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)
 	call dm_destroy(X,ierr)
 	call dm_destroy(Y,ierr)
 	call dm_destroy(Z,ierr)
 	call dm_destroy(U,ierr)


    if(myrank==0) print *, "==============Test dm_emult==============="
    A=dm_seqs(m,m)
    B=dm_eyes(m,m)
    C=A .em. B
    D=A .em. (dm_eyes(m,m))
    E=dm_seqs(m,m) .em. B
    F=dm_seqs(m,m) .em. dm_eyes(m,m) 
    G=A .em. A
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A.*B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=A.*dm_eyes(m,m*2)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_seqs(m,m).*B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=dm_seqs(m,m).*dm_eyes(m,m*2)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A.*A"
        call dm_view(G,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)


    if(myrank==0) print *, "==============Test dm_ediv================"
    A=dm_eyes(m,m)
    B=dm_seqs(m,m)
    C=A .ed. B
    D=A .ed. (dm_seqs(m,m))
    E=dm_eyes(m,m) .ed. B
    F=dm_eyes(m,m) .ed. dm_seqs(m,m) 
    G=A .ed. A
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A./B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=A./dm_seqs(m,m*2)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eyes(m,m)./B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=dm_eyes(m,m)./dm_seqs(m,m*2)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A./A"
        call dm_view(G,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)



    if(myrank==0) print *, "==============Test dm_rep================="
    A=dm_eyes(m,m)
    B=dm_rep(A,3,2) 
    C=dm_rep(dm_eyes(m,m),3,2) 
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_rep(A,3,2)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_rep(dm_eyes(m,m),3,2)"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_sum================="
    A=dm_seqs(m,n)
    B=dm_sum(A,1)
    C=dm_sum(A,2)
    D=dm_sum(dm_seqs(m,n),1)
    E=dm_sum(dm_seqs(m,n),2)
    if(debug) then
        if(myrank==0) print *, ">A=dm_seqs(m,n)"
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_sum(A,1)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_sum(A,2)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_sum(dm_seqs(m,n),2)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_sum(dm_seqs(m,n),2)"
        call dm_view(E,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)


    if(myrank==0) print *, "==============Test dm_axpy================"
    A=dm_seqs(m,m)	
    B=dm_eyes(m,m) 
    C=dm_eyes(m,m) 
    call dm_axpy(B,2.0,A,ierr) 
    call dm_axpy(C,2,dm_seqs(m,m),ierr) 
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B= dm_axpy(B,2.0,A)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_axpy(C,2,dm_seqs(m,m))"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_aypx================"
    A=dm_seqs(m,m)	
    B=dm_eyes(m,m) 
    C=dm_eyes(m,m) 
    call dm_aypx(B,2.0,A,ierr) 
    call dm_aypx(C,2,dm_seqs(m,m),ierr) 
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B= dm_axpy(B,2.0,A)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_axpy(C,2,dm_seqs(m,m))"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_trans==============="
    A=dm_seqs(m,n)	
    B=dm_trans(A)
    C=dm_trans(dm_seqs(m,n))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B= dm_trans(A)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_trans(dm_seqs(m,n))"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_xyt================="
    A=dm_seqs(m,m)
    B=dm_ones(m,m)
    C=dm_xyt(A,B)
    D=dm_xyt(A,dm_ones(m,m))
    E=dm_xyt(dm_seqs(m,m),B)
    F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))
    G=dm_xyt(A,A)
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_xyt(A,B)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_xyt(A,dm_ones(m,m))"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_xyt(dm_seqs(m,m),B)"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=dm_xyt(A,A)"
        call dm_view(G,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)


    if(myrank==0) print *, "==============Test dm_xty================="
    A=dm_seqs(m,m)
    B=dm_ones(m,m)
    C=dm_xty(A,B)
    D=dm_xty(A,dm_ones(m,m))
    E=dm_xty(dm_seqs(m,m),B)
    F=dm_xty(dm_seqs(m,m),dm_ones(m,m))
    G=dm_xty(A,A)
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_xyt(A,B)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_xyt(A,dm_ones(m,m))"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_xyt(dm_seqs(m,m),B)"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=dm_xyt(A,A)"
        call dm_view(G,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)

    if(myrank==0) print *, "==============Test dm_exp================="
    A=dm_seqs(m,m)
    B=dm_exp(A)
    C=dm_exp(dm_seqs(m,m))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_exp(A)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_exp(dm_seqs(m,m))"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_log================="
    A=dm_seqs(m,m)
    B=dm_log(A)
    C=dm_log(dm_seqs(m,m))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_log(A)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_log(dm_seqs(m,m))"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_sqrt================"
    A=dm_seqs(m,m)
    B=dm_sqrt(A)
    C=dm_sqrt(dm_seqs(m,m))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_sqrt(A)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_sqrt(dm_seqs(m,m))"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_squ================="
    A=dm_seqs(m,m)
    B=dm_squ(A)
    C=dm_squ(dm_seqs(m,m))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_squ(A)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_squ(dm_seqs(m,m))"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_cube================"
    A=dm_seqs(m,m)
    B=dm_cube(A)
    C=dm_cube(dm_seqs(m,m))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_cube(A)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_cube(dm_seqs(m,m))"
        call dm_view(C,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_solve==============="
    A=dm_seqs(m,m)
    B=dm_ones(m,1)
    C=dm_solve(A,B)
    D=A .inv. B
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_solve(A,B)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=A .inv. B"
        call dm_view(D,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)
  	call dm_destroy(D,ierr)


    if(myrank==0) print *, "==============Test dm_load================"
    filename="md001.00004"

    call dm_load(filename,A,ierr)
    if(debug) then
        if(myrank==0) print *, ">Load A from md001.00004="
        call dm_view(A,ierr)
 	endif
  	call dm_destroy(A,ierr)


    if(myrank==0) print *, "==============Test dm_setvalue============"
    A=dm_eyes(m,m)
    B=dm_eyes(m,m)
    C=dm_eyes(m,m)
    call dm_setvalue(A,1,1,5,ierr)
    call dm_setvalue(B,1,1,6.1,ierr)
    call dm_setvalue(C,1,1,real(7,8),ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C="
        call dm_view(C,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)

    if(myrank==0) print *, "==============Test dm_submatrix==========="
    A=dm_seqs(m,m)
    
    B=dm_zeros(2,1)
    call dm_setvalue(B,0,0,0,ierr)
    call dm_setvalue(B,1,0,2,ierr)
    
    C=dm_zeros(2,1)
    call dm_setvalue(C,0,0,0,ierr)
    call dm_setvalue(C,1,0,1,ierr)
    	
   	D=dm_submatrix(A,B,C)
    E=dm_submatrix(dm_seqs(m,m),B,C)
    F=dm_submatrix(A,dm_seqs(2,1),C)
    G=dm_submatrix(A,B,dm_seqs(2,1))
    H=dm_submatrix(dm_seqs(m,m),dm_seqs(2,1),C)
    X=dm_submatrix(dm_seqs(m,m),B,dm_seqs(2,1))
    Y=dm_submatrix(A,dm_seqs(2,1),dm_seqs(2,1))
    Z=dm_submatrix(dm_seqs(m,m),dm_seqs(2,1),dm_seqs(2,1))

    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C="
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_submatrix(A,B,C)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_submatrix(dm_seqs(m,m),B,C)"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=dm_submatrix(A,dm_seqs(2,1),C)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=dm_submatrix(A,B,dm_seqs(2,1))"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=dm_submatrix(dm_seqs(m,m),dm_seqs(2,1),C)"
        call dm_view(H,ierr)
        if(myrank==0) print *, ">X=dm_submatrix(dm_seqs(m,m),B,dm_seqs(2,1))"
        call dm_view(X,ierr)
        if(myrank==0) print *, ">Y=dm_submatrix(A,dm_seqs(2,1),dm_seqs(1,1))"
        call dm_view(Y,ierr)
        if(myrank==0) print *, ">Z=dm_submatrix(dm_seqs(m,m),dm_seqs(2,1),dm_seqs(2,1))"
        call dm_view(Z,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)
  	call dm_destroy(D,ierr)
  	call dm_destroy(E,ierr)
  	call dm_destroy(F,ierr)
  	call dm_destroy(G,ierr)
  	call dm_destroy(H,ierr)
  	call dm_destroy(X,ierr)
  	call dm_destroy(Y,ierr)
  	call dm_destroy(Z,ierr)


    if(myrank==0) print *, "==============Test dm_getcol=============="
    A=dm_seqs(m,n)
    B=dm_getcol(A,0)
    C=dm_getcol(A,1) 
    	
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_getcol(A,0)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_getcol(A,1)"
        call dm_view(C,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_getrow=============="
    A=dm_seqs(m,n)
    B=dm_getrow(A,0)
    C=dm_getrow(A,2) 
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_getrow(A,0)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_getrow(A,2)"
        call dm_view(C,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)


    if(myrank==0) print *, "==============Test dm_setvalues==========="
    A=dm_ones(2*m,2*m)
    idxm(1)=0	
    idxm(2)=2	
    idxn(1)=1	
    idxn(2)=2
    allocate(array(4))
    array(1)=9
    array(2)=8	
    array(3)=7	
    array(4)=6	
  	call dm_setvalues(A,2,idxm,2,idxn,array,ierr)	
    
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
 	endif
  	call dm_destroy(A,ierr)
    deallocate(array)

    if(myrank==0) print *, "==============Test dm_getvalues==========="
    A=dm_seqs(2*m,2*m)
    idxm(1)=A%ista	
    idxn(1)=0	
    idxn(2)=1
    allocate(array(2))
    array=0
  	call dm_getvalues(A,1,idxm,2,idxn,array,ierr)	
    
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">getvalues=",array
 	endif
  	call dm_destroy(A,ierr)
    deallocate(array)

    if(myrank==0) print *, "==============Test dm_norm================"
    A=dm_seqs(m,m)
   	a1=dm_norm_1(A)
    a2=dm_norm_2(A)
    a3=dm_norm_inf(A) 
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">NORM_1=",a1
        if(myrank==0) print *, ">NORM_2=",a2
        if(myrank==0) print *, ">NORM_INF=",a3
 	endif
  	call dm_destroy(A,ierr)

    if(myrank==0) print *, "==============Test dm_lt=================="
    A=dm_seqs(m,m)
    B=5*dm_ones(m,m)
    C=(A<B)
    D=(A<5)
    E=(A<5.0)
    F=(A<real(5.0,kind=8))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A<B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=A<5"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=A<5.0"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=A<real(5.0,kind=8)"
        call dm_view(F,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)
  	call dm_destroy(D,ierr)
  	call dm_destroy(E,ierr)
  	call dm_destroy(F,ierr)

    if(myrank==0) print *, "==============Test dm_le=================="
    A=dm_seqs(m,m)
    B=5*dm_ones(m,m)
    C=(A<=B)
    D=(A<=5)
    E=(A<=5.0)
    F=(A<=real(5.0,kind=8))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=(A<=B)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=(A<=5)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=(A<=5.0)"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=(A<=real(5.0,kind=8))"
        call dm_view(F,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)
  	call dm_destroy(D,ierr)
  	call dm_destroy(E,ierr)
  	call dm_destroy(F,ierr)

    if(myrank==0) print *, "==============Test dm_gt=================="
    A=dm_seqs(m,m)
    B=5*dm_ones(m,m)
    C=(A>B)
    D=(A>5)
    E=(A>5.0)
    F=(A>real(5.0,kind=8))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=(A>B)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=(A>5)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=(A>5.0)"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=(A>real(5.0,kind=8))"
        call dm_view(F,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)
  	call dm_destroy(D,ierr)
  	call dm_destroy(E,ierr)
  	call dm_destroy(F,ierr)

    if(myrank==0) print *, "==============Test dm_ge=================="
    A=dm_seqs(m,m)
    B=5*dm_ones(m,m)
    C=(A>=B)
    D=(A>=5)
    E=(A>=5.0)
    F=(A>=real(5.0,kind=8))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=(A>=B)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=(A>=5)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=(A>=5.0)"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=(A>=real(5.0,kind=8))"
        call dm_view(F,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)
  	call dm_destroy(D,ierr)
  	call dm_destroy(E,ierr)
  	call dm_destroy(F,ierr)

    if(myrank==0) print *, "==============Test dm_eq=================="
    A=dm_seqs(m,m)
    B=5*dm_ones(m,m)
    C=(A==B)
    D=(A==5)
    E=(A==5.0)
    F=(A==real(5.0,kind=8))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=(A==B)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=(A==5)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=(A==5.0)"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=(A==real(5.0,kind=8))"
        call dm_view(F,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)
  	call dm_destroy(D,ierr)
  	call dm_destroy(E,ierr)
  	call dm_destroy(F,ierr)


    if(myrank==0) print *, "==============Test dm_nq=================="
    A=dm_seqs(m,m)
    B=5*dm_ones(m,m)
    C=(A/=B)
    D=(A/=5)
    E=(A/=5.0)
    F=(A/=real(5.0,kind=8))
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B="
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=(A/=B)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=(A/=5)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=(A/=5.0)"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=(A/=real(5.0,kind=8))"
        call dm_view(F,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
  	call dm_destroy(C,ierr)
  	call dm_destroy(D,ierr)
  	call dm_destroy(E,ierr)
  	call dm_destroy(F,ierr)


    if(myrank==0) print *, "==============Test dm_cart2sph============"
    filename="md001.00004"
    call dm_load(filename,A,ierr)	
    call dm_cart2sph(A,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_cart2sph(A)"
        call dm_view(B,ierr)
 	endif
  	call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)


	call dm_finalize(ierr)
end program
