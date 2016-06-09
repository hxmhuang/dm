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
    real*8    		:: ep,alpha
	real*8			:: array(4)
	!real  			:: array(4)
    logical         :: debug 
    integer         :: ierr
    character(len=50):: filename
    debug=.false.

    ierr=dm_init()
    
    myrank=dm_comm_rank()
    
    mysize=dm_comm_size()
    
    m=dm_option_int('-m')
    n=dm_option_int('-n')
    ep=dm_option_real('-ep')

    debug=dm_option_bool('-debug')
	if(myrank==0) then 
       print *, "==============Input paramenters==========="
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

 	if(myrank==0) print *, "==============Test dm_seq================="
    A=dm_seqs(m,n)
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
 	endif
    ierr=dm_destroy(A)


    if(myrank==0) print *, "==============Test dm_eye================="
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


 	if(myrank==0) print *, "==============Test dm_copy================"
 	A=dm_eyes(m,m)
    B=A
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)


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
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=A+B"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_eyes(m,m)+dm_eyes(m,m)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_eyes(m,m)+B"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=A+dm_eyes(m,m)"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=A+A+A"
        ierr=dm_view(G)
        if(myrank==0) print *, ">H=B+G"
        ierr=dm_view(H)
        if(myrank==0) print *, ">X=A+2.0"
        ierr=dm_view(X)
        if(myrank==0) print *, ">Y=2+A"
        ierr=dm_view(Y)
        if(myrank==0) print *, ">Z=real(2,8)+A"
        ierr=dm_view(Z)
 	endif
     A=A+A	
     if(debug) then
         if(myrank==0) print *, ">A=A+A"
         ierr=dm_view(a)
 	 endif
    ierr=dm_destroy(A)
  	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
 	ierr=dm_destroy(F)
 	ierr=dm_destroy(G)
 	ierr=dm_destroy(H)
 	ierr=dm_destroy(X)
 	ierr=dm_destroy(Y)
 	ierr=dm_destroy(Z)


    if(myrank==0) print *, "==============Test dm_minus==============="
  	A=dm_ones(m,m)
    B=dm_eyes(m,m)
    C=A-B
    D=dm_eyes(m,m)-dm_eyes(m,m)
    E=dm_eyes(m,m)-B
    F=A-dm_eyes(m,m)
    G=A-A-A
    H=B-G
	X=A-2.0
	Y=2-A
	Z=real(2,8)-A
	U=(0-A)+2
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=A-B"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_eyes(m,m)-dm_eyes(m,m)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_eyes(m,m)-B"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=A-dm_eyes(m,m)"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=A-A-A"
        ierr=dm_view(G)
        if(myrank==0) print *, ">H=B-G"
        ierr=dm_view(H)
        if(myrank==0) print *, ">X=A-2.0"
        ierr=dm_view(X)
        if(myrank==0) print *, ">Y=2-A"
        ierr=dm_view(Y)
        if(myrank==0) print *, ">Z=real(2,8)-A"
        ierr=dm_view(Z)
        if(myrank==0) print *, ">U=(0-A)+2"
        ierr=dm_view(U)
 	endif
    A=A-A	
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
 	endif
    ierr=dm_destroy(A)
  	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
 	ierr=dm_destroy(F)
 	ierr=dm_destroy(G)
 	ierr=dm_destroy(H)
 	ierr=dm_destroy(X)
 	ierr=dm_destroy(Y)
 	ierr=dm_destroy(Z)
 	ierr=dm_destroy(U)


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
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=A .hjoin. B"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_eyes(m,m) .hjoin. dm_eyes(m,m)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_eyes(m,m) .hjoin. B"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=A .hjoin. dm_eyes(m,m)"
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


    if(myrank==0) print *, "==============Test dm_mult================"
    A=dm_eyes(m,m)
    B=dm_eyes(m,m*2)
    C=A*B
    D=A*(dm_eyes(m,m*2))
    E=dm_eyes(m,m)*B
    F=dm_eyes(m,m)*dm_eyes(m,m*2) 
    G=A*A
    X=A*2.0
    Y=2.0*A
    alpha=3.0
    Z=alpha*A
    U=A*alpha
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=A*B"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=A*dm_eyes(m,m*2)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_eyes(m,m)*B"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=dm_eyes(m,m)*dm_eyes(m,m*2)"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=A*A"
        ierr=dm_view(G)
        if(myrank==0) print *, ">X=2.0*A"
        ierr=dm_view(X)
        if(myrank==0) print *, ">Y=A*2.0"
        ierr=dm_view(Y)
        if(myrank==0) print *, ">Z=alpha*A"
        ierr=dm_view(Y)
        if(myrank==0) print *, ">U=A*alpha"
        ierr=dm_view(Y)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
    ierr=dm_destroy(F)
 	ierr=dm_destroy(G)
 	ierr=dm_destroy(X)
 	ierr=dm_destroy(Y)
 	ierr=dm_destroy(Z)
 	ierr=dm_destroy(U)


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


    if(myrank==0) print *, "==============Test dm_rep================="
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


    if(myrank==0) print *, "==============Test dm_sum================="
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


    if(myrank==0) print *, "==============Test dm_axpy================"
    A=dm_seqs(m,m)	
    B=dm_eyes(m,m) 
    C=dm_eyes(m,m) 
    alpha=1.0    
    ierr=dm_axpy(B,alpha,A) 
    ierr=dm_axpy(C,alpha,dm_seqs(m,m)) 
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B= dm_axpy(B,alpha,A)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_axpy(C,alpha,dm_seqs(m,m))"
        ierr=dm_view(C)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_aypx================"
    A=dm_seqs(m,m)	
    B=dm_eyes(m,m) 
    C=dm_eyes(m,m) 
    alpha=2.0    
    ierr=dm_aypx(B,alpha,A) 
    ierr=dm_aypx(C,alpha,dm_seqs(m,m)) 
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B= dm_axpy(B,alpha,A)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_axpy(C,alpha,dm_seqs(m,m))"
        ierr=dm_view(C)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_trans==============="
    A=dm_seqs(m,n)	
    B=dm_trans(A)
    C=dm_trans(dm_seqs(m,n))
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B= dm_trans(A)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_trans(dm_seqs(m,n))"
        ierr=dm_view(C)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)


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
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_xyt(A,B)"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_xyt(A,dm_ones(m,m))"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_xyt(dm_seqs(m,m),B)"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=dm_xyt(A,A)"
        ierr=dm_view(G)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
    ierr=dm_destroy(F)
 	ierr=dm_destroy(G)


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
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_xyt(A,B)"
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_xyt(A,dm_ones(m,m))"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_xyt(dm_seqs(m,m),B)"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=dm_xyt(A,A)"
        ierr=dm_view(G)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)
 	ierr=dm_destroy(D)
 	ierr=dm_destroy(E)
    ierr=dm_destroy(F)
 	ierr=dm_destroy(G)

    if(myrank==0) print *, "==============Test dm_exp================="
    A=dm_seqs(m,m)
    B=dm_exp(A)
    C=dm_exp(dm_seqs(m,m))
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B=dm_exp(A)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_exp(dm_seqs(m,m))"
        ierr=dm_view(C)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_log================="
    A=dm_seqs(m,m)
    B=dm_log(A)
    C=dm_log(dm_seqs(m,m))
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B=dm_log(A)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_log(dm_seqs(m,m))"
        ierr=dm_view(C)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_sqrt================"
    A=dm_seqs(m,m)
    B=dm_sqrt(A)
    C=dm_sqrt(dm_seqs(m,m))
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B=dm_sqrt(A)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_sqrt(dm_seqs(m,m))"
        ierr=dm_view(C)
 	endif
 	ierr=dm_destroy(A)
 	ierr=dm_destroy(B)
 	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_solve==============="
    A=dm_seqs(m,m)
    B=dm_ones(m,1)
    C=dm_solve(A,B)
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(b)
        if(myrank==0) print *, ">C=dm_solve(A,B)"
        ierr=dm_view(C)
 	endif
  	ierr=dm_destroy(A)
  	ierr=dm_destroy(B)
  	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_load================"
    filename="md001.00004"

    A=dm_load(filename)
    if(debug) then
        if(myrank==0) print *, ">Load A from md001.00004="
        ierr=dm_view(A)
 	endif
  	ierr=dm_destroy(A)


    if(myrank==0) print *, "==============Test dm_setvalue============"
    A=dm_eyes(m,m)
    B=dm_eyes(m,m)
    C=dm_eyes(m,m)
    ierr=dm_setvalue(A,1,1,5)
    ierr=dm_setvalue(B,1,1,6.1)
    ierr=dm_setvalue(C,1,1,real(7,8))
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C="
        ierr=dm_view(C)
 	endif
  	ierr=dm_destroy(A)
  	ierr=dm_destroy(B)
  	ierr=dm_destroy(C)

    if(myrank==0) print *, "==============Test dm_submatrix==========="
    A=dm_eyes(m,m)
    ierr=dm_setvalue(A,1,1,8)
    
    B=dm_ones(2,1)
    ierr=dm_setvalue(B,0,0,1)
    ierr=dm_setvalue(B,1,0,2)
    
    C=dm_ones(2,1)
    ierr=dm_setvalue(C,0,0,1)
    ierr=dm_setvalue(C,1,0,2)
    	
   	D=dm_submatrix(A,B,C)
    E=dm_submatrix(dm_eyes(m,m),B,C)
    F=dm_submatrix(A,dm_seqs(2,1),C)
    G=dm_submatrix(A,B,dm_seqs(2,1))
    H=dm_submatrix(dm_eyes(m,m),dm_seqs(2,1),C)
    X=dm_submatrix(dm_eyes(m,m),B,dm_seqs(2,1))
    Y=dm_submatrix(A,dm_seqs(2,1),dm_seqs(2,1))
    Z=dm_submatrix(dm_eyes(m,m),dm_seqs(2,1),dm_seqs(2,1))

    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B="
        ierr=dm_view(B)
        if(myrank==0) print *, ">C="
        ierr=dm_view(C)
        if(myrank==0) print *, ">D=dm_submatrix(A,B,C)"
        ierr=dm_view(D)
        if(myrank==0) print *, ">E=dm_submatrix(dm_eyes(m,m),B,C)"
        ierr=dm_view(E)
        if(myrank==0) print *, ">F=dm_submatrix(A,dm_seqs(2,1),C)"
        ierr=dm_view(F)
        if(myrank==0) print *, ">G=dm_submatrix(A,B,dm_seqs(2,1))"
        ierr=dm_view(G)
        if(myrank==0) print *, ">H=dm_submatrix(dm_eyes(m,m),dm_seqs(2,1),C)"
        ierr=dm_view(H)
        if(myrank==0) print *, ">X=dm_submatrix(dm_eyes(m,m),B,dm_seqs(2,1))"
        ierr=dm_view(X)
        if(myrank==0) print *, ">Y=dm_submatrix(A,dm_seqs(2,1),dm_seqs(1,1))"
        ierr=dm_view(Y)
        if(myrank==0) print *, ">Z=dm_submatrix(dm_eyes(m,m),dm_seqs(2,1),dm_seqs(2,1))"
        ierr=dm_view(Z)
 	endif
  	ierr=dm_destroy(A)
  	ierr=dm_destroy(B)
  	ierr=dm_destroy(C)
  	ierr=dm_destroy(D)
  	ierr=dm_destroy(E)
  	ierr=dm_destroy(F)
  	ierr=dm_destroy(G)
  	ierr=dm_destroy(H)
  	ierr=dm_destroy(X)
  	ierr=dm_destroy(Y)
  	ierr=dm_destroy(Z)


    if(myrank==0) print *, "==============Test dm_getcol=============="
    A=dm_seqs(m,n)
    B=dm_getcol(A,0)
    C=dm_getcol(A,1) 
    	
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B=dm_getcol(A,0)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_getcol(A,1)"
        ierr=dm_view(C)
 	endif
  	ierr=dm_destroy(A)
  	ierr=dm_destroy(B)
  	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_getrow=============="
    A=dm_seqs(m,n)
    B=dm_getrow(A,0)
    C=dm_getrow(A,2) 
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, ">B=dm_getrow(A,0)"
        ierr=dm_view(B)
        if(myrank==0) print *, ">C=dm_getrow(A,2)"
        ierr=dm_view(C)
 	endif
  	ierr=dm_destroy(A)
  	ierr=dm_destroy(B)
  	ierr=dm_destroy(C)


    if(myrank==0) print *, "==============Test dm_getsize=============="
    A=dm_eyes(m,m)
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
        if(myrank==0) print *, "nrow=",A%nrow,"ncol=",A%ncol
        if(myrank==0) print *, "rank0:ista=",A%ista,"iend=",A%iend
        if(myrank==1) print *, "rank1:ista=",A%ista,"iend=",A%iend
 	endif
  	ierr=dm_destroy(A)


   if(myrank==0) print *, "==============Test dm_setvalues============"
    A=dm_eyes(2*m,2*m)
    idxm(1)=0	
    idxm(2)=2	
    idxn(1)=1	
    idxn(2)=2
    array(1)=9
    array(2)=8	
    array(3)=7	
    array(4)=6	
  	ierr=dm_setvalues(A,2,idxm,2,idxn,array)	
    
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
 	endif
  	ierr=dm_destroy(A)


	ierr=dm_finalize()
end program
