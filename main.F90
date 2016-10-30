program main 

    use dm 
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
    type(Matrix)    		:: A,B,C,D,E,F,G,H 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: myrank, mysize 
    integer         		:: m,n,k
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug 
    integer         		:: ierr
    character(len=50)		:: filename
	real(kind=8),allocatable :: array(:)
    debug=.false.

    call dm_init1(ierr)
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
	call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
    
    call dm_init2(m,n,k,ierr)
    if(myrank==0) then 
       print *, "==============Input paramenters==========="
        print *, "m=",m,",n=",n,",k=",k,",ep=",ep,",debug=",debug
     endif 
    
    
 	if(myrank==0) print *, "==============Test dm_zeros==============="
    A=dm_zeros(m,n,k)
    B=dm_zeros(m,n,k,.true.)
    C=dm_zeros(m,n,k,.false.)
    D=dm_zeros(m,m,m)
    if(debug) then
        if(myrank==0) print *, ">A=dm_zeros(m,n,k)"
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_zeros(m,n,k,.true.)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_zeros(m,n,k,.false.)"
        if(myrank==0) call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_zeros(m,m,m)"
        call dm_view(D,ierr)
 	endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)


 	if(myrank==0) print *, "==============Test dm_ones================"
    A=dm_ones(m,n,k)
    B=dm_ones(m,n,k,.true.)
    C=dm_ones(m,n,k,.false.)
    D=dm_ones(m,m,m)
    E=dm_ones(n,n,m)
    if(debug) then
        if(myrank==0) print *, ">A=dm_ones(m,n,k)"
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_ones(m,n,k,.true.)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_ones(m,n,k,.false.)"
        if(myrank==0) call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_ones(m,m,m)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_ones(n,n,m)"
        call dm_view(E,ierr)
 	endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)


    if(myrank==0) print *, "==============Test dm_eye================="
    A=dm_eye(m,n,k)	
    B=dm_eye(m,m*2,1)	
    C=dm_eye(2*m,m,1)	
    D=dm_eye(m,m*2,1,.true.)	
    E=dm_eye(2*m,m,1,.false.)	
    F=dm_eye(m,m*2,2)	
    G=dm_eye(2*m,m,3)
    H=dm_eye(m,m,m)	
    if(debug) then
        if(myrank==0) print *, ">A=dm_eye(m,n,k)"
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_eye(m,2m,1)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=dm_eye(2*m,m,1)"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_eye(m,m*2,1,.true.)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eye(2*m,m,1,.false.)"
        if(myrank==0) call dm_view(E,ierr)
        if(myrank==0) print *, ">F=dm_eye(m,2*m,2)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=dm_eye(2*m,m,3)"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=dm_eye(m,m,m)"
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


 	if(myrank==0) print *, "==============Test dm_copy================"
 	A=dm_eye(m,n,k)
    B=A
    if(debug) then
        if(myrank==0) print *, ">A=dm_eye(m,n,k)"
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=A"
        call dm_view(B,ierr)
 	endif
 	call dm_destroy(A,ierr)
 	call dm_destroy(B,ierr)


    if(myrank==0) print *, "==============Test dm_add================="
  	A=dm_eye(m,n,k)
    B=A
    C=A+B
    D=dm_eye(m,n,k)+dm_eye(m,n,k)
    E=dm_eye(m,n,k)+B
    F=A+dm_eye(m,n,k)
    G=A+A+A
    H=B+G
    X=A+2.0
    Y=2+A
    Z=real(2,8)+A
    U=dm_eye(m,n,k,.true.)+dm_eye(m,n,k)
    V=dm_eye(m,n,k,.false.)+dm_eye(m,n,k,.false.)
    W=2+dm_eye(m,n,k,.false.)
    if(debug) then
        if(myrank==0) print *, ">A=dm_eye(m,n,k)"
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=A"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A+B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_eye(m,n,k)+dm_eye(m,n,k)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eye(m,n,k)+B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=A+dm_eye(m,n,k)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A+A+A"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=B+G"
        call dm_view(H,ierr)
        if(myrank==0) print *, ">X=A+2.0"
        call dm_view(X,ierr)
        if(myrank==0) print *, ">Y=2+A"
        call dm_view(Y,ierr)
        if(myrank==0) print *, ">Z=real(2,kind=8)+A"
        call dm_view(Z,ierr)
        if(myrank==0) print *, ">U=dm_eye(m,n,k,.true.)+dm_eye(m,n,k)"
        call dm_view(U,ierr)
        if(myrank==0) print *, ">V=dm_eye(m,n,k,.false.)+dm_eye(m,n,k,.false.)"
        if(myrank==0) call dm_view(V,ierr)
        if(myrank==0) print *, ">W=2+dm_eye(m,n,k,.false.)"
        if(myrank==0) call dm_view(W,ierr)
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
 	call dm_destroy(U,ierr)
 	call dm_destroy(V,ierr)
 	call dm_destroy(W,ierr)


    if(myrank==0) print *, "==============Test dm_minus==============="
  	A=dm_ones(m,n,k)
    B=dm_eye(m,n,k)
    C=A-B
    D=dm_eye(m,n,k)-dm_eye(m,n,k)
    E=dm_eye(m,n,k)-B
    F=A-dm_eye(m,n,k)
    G=A-A-A
    H=B-G
    X=A-2.0
    Y=2-A
    Z=real(2,8)-A
    U=(0-A)+2
    V=dm_eye(m,n,k,.false.)-dm_eye(m,n,k,.false.)
    W=2-dm_eye(m,n,k,.false.)
    if(debug) then
        if(myrank==0) print *, ">A=dm_ones(m,n,k)"
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm_eye(m,n,k)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A-B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_eye(m,n,k)-dm_eye(m,n,k)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eye(m,n,k)-B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=A-dm_eye(m,n,k)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A-A-A"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=B-G"
        call dm_view(H,ierr)
        if(myrank==0) print *, ">X=A-2.0"
        call dm_view(X,ierr)
        if(myrank==0) print *, ">Y=2-A"
        call dm_view(Y,ierr)
        if(myrank==0) print *, ">Z=real(2,8)-A"
        call dm_view(Z,ierr)
        if(myrank==0) print *, ">U=(0-A)+2"
        call dm_view(U,ierr)
        if(myrank==0) print *, ">V=dm_eye(m,n,k,.false.)-dm_eye(m,n,k,.false.)"
        if(myrank==0) call dm_view(V,ierr)
        if(myrank==0) print *, ">W=2-dm_eye(m,n,k,.false.)"
        if(myrank==0) call dm_view(W,ierr)
 	endif
    A=-A	
    if(debug) then
        if(myrank==0) print *, ">A=-A"
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
 	call dm_destroy(V,ierr)
 	call dm_destroy(W,ierr)


    if(myrank==0) print *, "==============Test dm_xjoin==============="
  	A=dm_eye(m,n,k)
    B=dm_eye(m,n,k)
    C=A .xj. B
    D=dm_eye(m,n,k) .xj. dm_eye(m,n,k)
    E=dm_eye(m,n,k) .xj. B
    F=A .xj. dm_eye(m,n,k)
    G=A .xj. A .xj. A
    H=B .xj. G
    U=dm_eye(m,n,k,.false.) .xj. dm_eye(m,n,k,.false.)
    if(debug) then
        if(myrank==0) print *, ">A=dm(m,n,k)"
        call dm_view(A,ierr)
        if(myrank==0) print *, ">B=dm(m,n,k)"
        call dm_view(B,ierr)
        if(myrank==0) print *, ">C=A .xj. B"
        call dm_view(C,ierr)
        if(myrank==0) print *, ">D=dm_eye(m,n,k) .xj. dm_eye(m,n,k)"
        call dm_view(D,ierr)
        if(myrank==0) print *, ">E=dm_eye(m,n,k) .xj. B"
        call dm_view(E,ierr)
        if(myrank==0) print *, ">F=A .xj. dm_eye(m,n,k)"
        call dm_view(F,ierr)
        if(myrank==0) print *, ">G=A .xj. A .xj. A"
        call dm_view(G,ierr)
        if(myrank==0) print *, ">H=B .xj. G"
        call dm_view(H,ierr)
        if(myrank==0) print *, ">U=dm_eye(m,n,k,.false.) .xj. dm_eye(m,n,k,.false.)"
        if(myrank==0) call dm_view(U,ierr)
 	endif
    call dm_destroy(A,ierr)
  	call dm_destroy(B,ierr)
 	call dm_destroy(C,ierr)
 	call dm_destroy(D,ierr)
 	call dm_destroy(E,ierr)
 	call dm_destroy(F,ierr)
 	call dm_destroy(G,ierr)
 	call dm_destroy(H,ierr)
 	call dm_destroy(U,ierr)


!   if(myrank==0) print *, "==============Test dm_vjoin==============="
! 	A=dm_eye(m,n,k)
!   B=dm_eye(m,n,k)
!   C=A .vj. B
!   D=dm_eye(m,n,k) .vj. dm_eye(m,n,k)
!   E=dm_eye(m,n,k) .vj. B
!   F=A .vj. dm_eye(m,n,k)
!   G=A .vj. A .vj. A
!   H=B .vj. G
!   U=dm_eye(m,n,k,.false.) .vj. dm_eye(m,n,k,.false.)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=A .vj. B"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_eye(m,n,k) .vj. dm_eye(m,n,k)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_eye(m,n,k) .vj. B"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=A .vj. dm_eye(m,n,k)"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=A .vj. A .vj. A"
!       call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=B .vj. G"
!       call dm_view(H,ierr)
!       if(myrank==0) print *, ">U=dm_eye(m,n,k,.false.) .vj. dm_eye(m,n,k,.false.)"
!       if(myrank==0) call dm_view(U,ierr)
!	endif
!   call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)
!	call dm_destroy(E,ierr)
!	call dm_destroy(F,ierr)
!	call dm_destroy(G,ierr)
!	call dm_destroy(H,ierr)
!	call dm_destroy(U,ierr)


!   if(myrank==0) print *, "==============Test dm_mult================"
!   A=dm_eye(m,n,k)
!   B=dm_eye(m,m*2)
!   C=A*B
!   D=A*(dm_eye(m,m*2))
!   E=dm_eye(m,n,k)*B
!   F=dm_eye(m,n,k)*dm_eye(m,m*2) 
!   G=A*A
!   X=A*2.0
!   Y=2*A
!   alpha=3.0
!   Z=alpha*A
!   U=A*alpha
!   V=dm_eye(m,n,k,.false.)*dm_eye(m,n,k,.false.)
!   W=3*dm_eye(m,n,k,.false.)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=A*B"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=A*dm_eye(m,m*2)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_eye(m,n,k)*B"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_eye(m,n,k)*dm_eye(m,m*2)"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=A*A"
!       call dm_view(G,ierr)
!       if(myrank==0) print *, ">X=A*2.0"
!       call dm_view(X,ierr)
!       if(myrank==0) print *, ">Y=2*A"
!       call dm_view(Y,ierr)
!       if(myrank==0) print *, ">Z=alpha*A"
!       call dm_view(Y,ierr)
!       if(myrank==0) print *, ">U=A*alpha"
!       call dm_view(Y,ierr)
!       if(myrank==0) print *, ">V=dm_eye(m,n,k,.false.)*dm_eye(m,n,k,.false.)"
!       if(myrank==0) call dm_view(V,ierr)
!       if(myrank==0) print *, ">W=3*dm_eye(m,n,k,.false.)"
!       if(myrank==0) call dm_view(W,ierr)

!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)
!	call dm_destroy(E,ierr)
!   call dm_destroy(F,ierr)
!	call dm_destroy(G,ierr)
!	call dm_destroy(X,ierr)
!	call dm_destroy(Y,ierr)
!	call dm_destroy(Z,ierr)
!	call dm_destroy(U,ierr)
!	call dm_destroy(V,ierr)
!	call dm_destroy(W,ierr)


!   if(myrank==0) print *, "==============Test dm_emult==============="
!   A=dm_seqs(m,m)
!   B=dm_eye(m,n,k)
!   C=A .em. B
!   D=A .em. (dm_eye(m,n,k))
!   E=dm_seqs(m,m) .em. B
!   F=dm_seqs(m,m) .em. dm_eye(m,n,k) 
!   G=A .em. A
!   H=dm_eye(m,n,k,.false.) .em. dm_seqs(m,m,.false.)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=A.*B"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=A.*dm_eye(m,m*2)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_seqs(m,m).*B"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_seqs(m,m).*dm_eye(m,m*2)"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=A.*A"
!       call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=dm_eye(m,n,k,.false.) .em. dm_seqs(m,m,.false.)"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)
!	call dm_destroy(E,ierr)
!   call dm_destroy(F,ierr)
!	call dm_destroy(G,ierr)
!	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_ediv================"
!   A=dm_eye(m,n,k)
!   B=dm_seqs(m,m)
!   C=A .ed. B
!   D=A .ed. (dm_seqs(m,m))
!   E=dm_eye(m,n,k) .ed. B
!   F=dm_eye(m,n,k) .ed. dm_seqs(m,m) 
!   G=A .ed. A
!   H=dm_eye(m,n,k,.false.) .ed. dm_eye(m,n,k,.false.)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=A./B"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=A./dm_seqs(m,m*2)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_eye(m,n,k)./B"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_eye(m,n,k)./dm_seqs(m,m*2)"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=A./A"
!       call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=dm_eye(m,n,k,.false.) .ed. dm_eye(m,n,k,.false.)"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)
!	call dm_destroy(E,ierr)
!   call dm_destroy(F,ierr)
!	call dm_destroy(G,ierr)
!	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_rep================="
!   A=dm_eye(m,n,k)
!   B=dm_rep(A,3,2) 
!   C=dm_rep(dm_eye(m,n,k),3,2) 
!   D=dm_rep(dm_eye(m,n,k,.false.),3,2) 
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_rep(A,3,2)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_rep(dm_eye(m,n,k),3,2)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_rep(dm_eye(m,n,k,.false.),3,2)"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_sum================="
!   A=dm_seqs(m,n)
!   B=dm_sum(A,1)
!   C=dm_sum(A,2)
!   D=dm_sum(dm_seqs(m,n),1)
!   E=dm_sum(dm_seqs(m,n),2)
!   F=dm_sum(dm_seqs(m,n,.false.),2)
!   if(debug) then
!       if(myrank==0) print *, ">A=dm_seqs(m,n)"
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_sum(A,1)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_sum(A,2)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_sum(dm_seqs(m,n),2)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_sum(dm_seqs(m,n),2)"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_sum(dm_seqs(m,n.false.),2)"
!       if(myrank==0) call dm_view(F,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)
!	call dm_destroy(E,ierr)
!	call dm_destroy(F,ierr)


!   if(myrank==0) print *, "==============Test dm_axpy================"
!   A=dm_seqs(m,m)	
!   B=dm_eye(m,n,k) 
!   C=dm_eye(m,n,k) 
!   D=dm_eye(m,n,k,.false.) 

!   call dm_axpy(B,2.0,A,ierr) 
!   call dm_axpy(C,2,dm_seqs(m,m),ierr) 
!   call dm_axpy(D,2,dm_seqs(m,m,.false.),ierr) 
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B= dm_axpy(B,2.0,A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_axpy(C,2,dm_seqs(m,m))"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_axpy(D,2,dm_seqs(m,m,.false.))"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_aypx================"
!   A=dm_seqs(m,m)	
!   B=dm_eye(m,n,k) 
!   C=dm_eye(m,n,k) 
!   D=dm_eye(m,n,k,.false.) 
!   call dm_aypx(B,2.0,A,ierr) 
!   call dm_aypx(C,2,dm_seqs(m,m),ierr) 
!   call dm_aypx(D,2,dm_seqs(m,m,.false.),ierr) 
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B= dm_axpy(B,2.0,A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_axpy(C,2,dm_seqs(m,m))"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_axpy(D,2,dm_seqs(m,m,.false.))"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)


!   if(myrank==0) print *, "==============Test dm_trans==============="
!   A=dm_seqs(m,n)	
!   B=dm_trans(A)
!   C=dm_trans(dm_seqs(m,n))
!   D=dm_trans(dm_seqs(m,n,.false.))
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B= dm_trans(A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_trans(dm_seqs(m,n))"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_trans(dm_seqs(m,m,.false.))"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_xyt================="
!   A=dm_seqs(m,m)
!   B=dm_ones(m,m)
!   C=dm_xyt(A,B)
!   D=dm_xyt(A,dm_ones(m,m))
!   E=dm_xyt(dm_seqs(m,m),B)
!   F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))
!   G=dm_xyt(A,A)
!   H=dm_xyt(dm_seqs(m,m,.false.),dm_ones(m,m,.false.))
!   if(debug) then
!       if(myrank==0) print *, ">A=dm_seqs(m,m)"
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_ones(m,m)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_xyt(A,B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_xyt(A,dm_ones(m,m))"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_xyt(dm_seqs(m,m),B)"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=dm_xyt(A,A)"
!       call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=dm_xyt(dm_seqs(m,m,.false.),dm_ones(m,m,.false.))"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)
!	call dm_destroy(E,ierr)
!   call dm_destroy(F,ierr)
!	call dm_destroy(G,ierr)
!	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_xty================="
!   A=dm_seqs(m,m)
!   B=dm_ones(m,m)
!   C=dm_xty(A,B)
!   D=dm_xty(A,dm_ones(m,m))
!   E=dm_xty(dm_seqs(m,m),B)
!   F=dm_xty(dm_seqs(m,m),dm_ones(m,m))
!   G=dm_xty(A,A)
!   H=dm_xty(dm_seqs(m,m,.false.),dm_ones(m,m,.false.))
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_xyt(A,B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_xyt(A,dm_ones(m,m))"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_xyt(dm_seqs(m,m),B)"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=dm_xyt(A,A)"
!       call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=dm_xty(dm_seqs(m,m,.false.),dm_ones(m,m,.false.))"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)
!	call dm_destroy(E,ierr)
!   call dm_destroy(F,ierr)
!	call dm_destroy(G,ierr)
!	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_exp================="
!   A=dm_seqs(m,m)
!   B=dm_exp(A)
!   C=dm_exp(dm_seqs(m,m))
!   D=dm_exp(dm_seqs(m,m,.false.))
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_exp(A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_exp(dm_seqs(m,m))"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_exp(dm_seqs(m,m,.false.))"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_log================="
!   A=dm_seqs(m,m)
!   B=dm_log(A)
!   C=dm_log(dm_seqs(m,m))
!   D=dm_log(dm_seqs(m,m,.false.))
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_log(A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_log(dm_seqs(m,m))"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_log(dm_seqs(m,m,.false.))"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_sqrt================"
!   A=dm_seqs(m,m)
!   B=dm_sqrt(A)
!   C=dm_sqrt(dm_seqs(m,m))
!   D=dm_sqrt(dm_seqs(m,m,.false.))
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_sqrt(A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_sqrt(dm_seqs(m,m))"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_sqrt(dm_seqs(m,m,.false.))"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_squ================="
!   A=dm_seqs(m,m)
!   B=dm_squ(A)
!   C=dm_squ(dm_seqs(m,m))
!   D=dm_squ(dm_seqs(m,m,.false.))
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_squ(A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_squ(dm_seqs(m,m))"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_squ(dm_seqs(m,m,.false.))"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_cube================"
!   A=dm_seqs(m,m)
!   B=dm_cube(A)
!   C=dm_cube(dm_seqs(m,m))
!   D=dm_cube(dm_seqs(m,m,.false.))
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_cube(A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_cube(dm_seqs(m,m))"
!       if(myrank==0) print *, ">D=dm_cube(dm_seqs(m,m,.false.))"
!       if(myrank==0) call dm_view(D,ierr)
!       call dm_view(C,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
!	call dm_destroy(C,ierr)
!	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_solve==============="
!   A=dm_seqs(m,m)
!   B=dm_ones(m,1)
!   C=dm_solve(A,B)
!   D=A .inv. B
!   !E=1+dm_seqs(m,m,.false.)
!   E=dm_ones(m,m,.false.)
!   F=3*dm_ones(m,1,.false.)
!   G=dm_solve(E,F)
!   H=E .inv. F
!   if(debug) then
!       if(myrank==0) print *, ">A=dm_seqs(m,m)"
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_ones(m,1)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_solve(A,B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=A .inv. B"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_seqs(m,m,.false.)"
!       if(myrank==0) call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=3*dm_ones(m,1,.false.)"
!       if(myrank==0) call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=dm_solve(E,F)"
!       if(myrank==0) call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=E .inv. F"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_load================"
!   filename="md001.00004"

!   call dm_load(filename,.true.,A,ierr)
!   call dm_load(filename,.false.,B,ierr)
!   if(debug) then
!       if(myrank==0) print *, ">Load A from md001.00004="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">Load local B from md001.00004="
!       if(myrank==0) call dm_view(B,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)


!   if(myrank==0) print *, "==============Test dm_setvalue============"
!   A=dm_eye(m,n,k)
!   B=dm_eye(m,n,k)
!   C=dm_eye(m,n,k)
!   D=dm_eye(m,n,k,.false.)
!   call dm_setvalue(A,1,1,5,ierr)
!   call dm_setvalue(B,1,1,6.1,ierr)
!   call dm_setvalue(C,1,1,real(7,kind=8),ierr)
!   call dm_setvalue(D,1,2,8,ierr)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C="
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D="
!       if(myrank==0) call dm_view(D,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)

!   if(myrank==0) print *, "==============Test dm_getsub=============="
!   A=dm_seqs(m,m)
!  	B=dm_getsub(A,(/0,1,2/),(/0,1/))
!  	C=dm_getsub(A,(/0,2,1/),(/1,0/))
!   D=dm_seqs(m,m,.false.)
!  	E=dm_getsub(D,(/0,1,2/),(/0,1/))
!  	F=dm_getsub(D,(/0,2,1/),(/1,0/))
!   if(debug) then
!       if(myrank==0) print *, ">A=dm_seqs(m,m)"
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_getsub(A,(/0,1,2/),(/0,1/))"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_getsub(A,(/0,2,1/),(/1,0/))"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_seqs(m,m,.false.)"
!       if(myrank==0) call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_getsub(C,(/0,1,2/),(/0,1/))"
!       if(myrank==0) call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_getsub(C,(/0,1,2/),(/0,1/))"
!       if(myrank==0) call dm_view(F,ierr)
!	endif
!	call dm_destroy(A,ierr)
!	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)


!   if(myrank==0) print *, "==============Test dm_getcol=============="
!   A=dm_seqs(m,m)
!   B=dm_getcol(A,0)
!   C=dm_getcol(A,1) 
!   D=dm_seqs(m,m,.false.)
!   E=dm_getcol(D,0)
!   F=dm_getcol(D,1) 
!   if(debug) then
!       if(myrank==0) print *, ">A=dm_seqs(m,m)"
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_getcol(A,0)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_getcol(A,1)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_seqs(m,m,.false.)"
!       if(myrank==0) call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_getcol(D,0)"
!       if(myrank==0) call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_getcol(D,1)"
!       if(myrank==0) call dm_view(F,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)


!   if(myrank==0) print *, "==============Test dm_getrow=============="
!   A=dm_seqs(m,m)
!   B=dm_getrow(A,0)
!   C=dm_getrow(A,2) 
!   D=dm_seqs(m,m,.false.)
!   E=dm_getrow(D,0)
!   F=dm_getrow(D,1) 
!   if(debug) then
!       if(myrank==0) print *, ">A=dm_seqs(m,m)"
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_getrow(A,0)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_getrow(A,2)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_seqs(m,m,.false.)"
!       if(myrank==0) call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_getrow(D,0)"
!       if(myrank==0) call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_getrow(D,1)"
!       if(myrank==0) call dm_view(F,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)


!   if(myrank==0) print *, "==============Test dm_setvalues==========="
!   A=dm_ones(2*m,2*m)
!   B=dm_ones(2*m,2*m,.false.)
! 	call dm_setvalues(A,(/0,2/),(/1,2/),(/9,8,7,6/),ierr)	
! 	call dm_setvalues(B,(/0,2/),(/1,2/),(/9,8,7,6/),ierr)	
!   
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       if(myrank==0) call dm_view(B,ierr)
!	endif
! 	call dm_destroy(A,ierr)

!   if(myrank==0) print *, "==============Test dm_getvalues==========="
!   A=dm_seqs(2*m,2*m)
!   B=dm_seqs(2*m,2*m,.false.)
!   allocate(array(2))
!   array=0
! 	call dm_getvalues(A,(/A%ista/),(/1,2/),array,ierr)	
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">getvalues=dm_getvalues(A,idxm,idxn,array,ierr) ",array
!	endif
! 	call dm_getvalues(B,(/A%ista/),(/1,2/),array,ierr)	
!   if(debug) then
!       if(myrank==0) print *, ">B="
!       if(myrank==0) call dm_view(B,ierr)
!       if(myrank==0) print *, ">getvalues=dm_getvalues(B,idxm,idxn,array,ierr)",array
!	endif
!   call dm_destroy(B,ierr)
!	deallocate(array)

!   if(myrank==0) print *, "==============Test dm_norm================"
!   A=dm_seqs(m,m)
!  	a1=dm_norm_1(A)
!   a2=dm_norm_2(A)
!   a3=dm_norm_inf(A) 
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">NORM_1=",a1
!       if(myrank==0) print *, ">NORM_2=",a2
!       if(myrank==0) print *, ">NORM_INF=",a3
!	endif
!   B=dm_seqs(m,m,.false.)
!  	a1=dm_norm_1(B)
!   a2=dm_norm_2(B)
!   a3=dm_norm_inf(B) 
!   if(debug) then
!       if(myrank==0) print *, ">B="
!       if(myrank==0) call dm_view(B,ierr)
!       if(myrank==0) print *, ">NORM_1=",a1
!       if(myrank==0) print *, ">NORM_2=",a2
!       if(myrank==0) print *, ">NORM_INF=",a3
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)

!   if(myrank==0) print *, "==============Test dm_lt=================="
!   A=dm_seqs(m,m)
!   B=5*dm_ones(m,m)
!   C=(A<B)
!   D=(A<5)
!   E=(A<5.0)
!   F=(A<real(5.0,kind=8))
!   G=(dm_seqs(m,m,.false.) < dm_constants(m,m,5,.false.))
!   H=(dm_seqs(m,m,.false.) < 5)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=A<B"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=A<5"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=A<5.0"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=A<real(5.0,kind=8)"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=(dm_seqs(m,m,.false.) < dm_constants(m,m,5,.false.))"
!       if(myrank==0) call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=(dm_seqs(m,m,.false.) < 5)"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_le=================="
!   A=dm_seqs(m,m)
!   B=5*dm_ones(m,m)
!   C=(A<=B)
!   D=(A<=5)
!   E=(A<=5.0)
!   F=(A<=real(5.0,kind=8))
!   G=(dm_seqs(m,m,.false.) <= dm_constants(m,m,5,.false.))
!   H=(dm_seqs(m,m,.false.) <= 5)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=(A<=B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=(A<=5)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=(A<=5.0)"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=(A<=real(5.0,kind=8))"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=(dm_seqs(m,m,.false.) <= dm_constants(m,m,5,.false.))"
!       if(myrank==0) call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=(dm_seqs(m,m,.false.) <= 5)"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)

!   if(myrank==0) print *, "==============Test dm_gt=================="
!   A=dm_seqs(m,m)
!   B=5*dm_ones(m,m)
!   C=(A>B)
!   D=(A>5)
!   E=(A>5.0)
!   F=(A>real(5.0,kind=8))
!   G=(dm_seqs(m,m,.false.) > dm_constants(m,m,5,.false.))
!   H=(dm_seqs(m,m,.false.) > 5)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=(A>B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=(A>5)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=(A>5.0)"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=(A>real(5.0,kind=8))"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=(dm_seqs(m,m,.false.) > dm_constants(m,m,5,.false.))"
!       if(myrank==0) call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=(dm_seqs(m,m,.false.) > 5)"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)

!   if(myrank==0) print *, "==============Test dm_ge=================="
!   A=dm_seqs(m,m)
!   B=5*dm_ones(m,m)
!   C=(A>=B)
!   D=(A>=5)
!   E=(A>=5.0)
!   F=(A>=real(5.0,kind=8))
!   G=(dm_seqs(m,m,.false.) >= dm_constants(m,m,5,.false.))
!   H=(dm_seqs(m,m,.false.) >= 5)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=(A>=B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=(A>=5)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=(A>=5.0)"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=(A>=real(5.0,kind=8))"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=(dm_seqs(m,m,.false.) >= dm_constants(m,m,5,.false.))"
!       if(myrank==0) call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=(dm_seqs(m,m,.false.) >= 5)"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)

!   if(myrank==0) print *, "==============Test dm_eq=================="
!   A=dm_seqs(m,m)
!   B=5*dm_ones(m,m)
!   C=(A==B)
!   D=(A==5)
!   E=(A==5.0)
!   F=(A==real(5.0,kind=8))
!   G=(dm_seqs(m,m,.false.) == dm_constants(m,m,5,.false.))
!   H=(dm_seqs(m,m,.false.) == 5)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=(A==B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=(A==5)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=(A==5.0)"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=(A==real(5.0,kind=8))"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=(dm_seqs(m,m,.false.) == dm_constants(m,m,5,.false.))"
!       if(myrank==0) call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=(dm_seqs(m,m,.false.) == 5)"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_nq=================="
!   A=dm_seqs(m,m)
!   B=5*dm_ones(m,m)
!   C=(A/=B)
!   D=(A/=5)
!   E=(A/=5.0)
!   F=(A/=real(5.0,kind=8))
!   G=(dm_seqs(m,m,.false.) /= dm_constants(m,m,5,.false.))
!   H=(dm_seqs(m,m,.false.) /= 5)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=(A/=B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=(A/=5)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=(A/=5.0)"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=(A/=real(5.0,kind=8))"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=(dm_seqs(m,m,.false.) /= dm_constants(m,m,5,.false.))"
!       if(myrank==0) call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=(dm_seqs(m,m,.false.) /= 5)"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_sparse=============="
!   A=dm_m2n(1,3) .vj. dm_m2n(5,6)
!   B=dm_m2n(0,2) .vj. dm_m2n(5,6)
!   C=dm_m2n(100,104)
!   D=dm_sparse(A,B,C,8,8)	
!   E=dm_m2n(1,3,.false.) .vj. dm_m2n(5,6,.false.)
!   F=dm_m2n(0,2,.false.) .vj. dm_m2n(5,6,.false.)
!   G=dm_m2n(100,104,.false.)
!   H=dm_sparse(E,F,G,8,8)	
!   if(debug) then
!       if(myrank==0) print *, ">Ind_i:A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">Ind_j:B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_seqs(m,1)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_sparse"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">Ind_i:E="
!       if(myrank==0) call dm_view(E,ierr)
!       if(myrank==0) print *, ">Ind_j:F="
!       if(myrank==0) call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=dm_seqs(m,1)"
!       if(myrank==0) call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=dm_sparse"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_cart2sph============"
!   filename="md001.00004"
!   call dm_load(filename,.true.,A,ierr)	
!   call dm_load(filename,.false.,C,ierr)	
!   call dm_cart2sph(A,B,ierr)
!   call dm_cart2sph(C,D,ierr)
!   if(debug) then
!       if(myrank==0) print *, ">A=dm_load(filename,A,.true.ierr) "
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_cart2sph(A)"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_load(filename,A,.false.ierr) "
!       if(myrank==0) call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_cart2sph(A)"
!       if(myrank==0) call dm_view(D,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)


!   if(myrank==0) print *, "==============Test dm_setdiag============="
!   A=dm_seqs(m,m)
!   B=A
!   C=A
!   D=A
!   E=dm_seqs(m,m,.false.)
!   call dm_setdiag(B,real(2.0,kind=8),ierr)
!   call dm_setdiag(C,1.5,ierr)
!   call dm_setdiag(D,1,ierr)
!   call dm_setdiag(E,1,ierr)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B=dm_set_diag(A,real(2.0,kind=8))"
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_set_diag(A,1.5)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_set_diag(A,1)"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_set_diag(E,1)"
!       if(myrank==0) call dm_view(E,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)


!   if(myrank==0) print *, "==============Test dm_setcol=============="
!   A=dm_seqs(m,m)
!   !B=dm_m2n(10,10+m-1)
!   B=dm_veyezero(1,m-1)
!   C=A
!   D=A
!   E=A
!   F=A
!   G=A
!   H=dm_seqs(m,m,.false.)
!   call dm_setcol(C,0,B,ierr)
!	call dm_setcol(D,m-1,B,ierr)
!	call dm_setcol(E,m-1,dm_m2n(10,10+m-1),ierr)
!	call dm_setcol(F,1,dm_zeros(m,1),ierr)
!	call dm_setcol(G,1,dm_ones(m,1),ierr)
!	call dm_setcol(H,1,dm_ones(m,1,.false.),ierr)

!   if(debug) then
!       if(myrank==0) print *, ">A="
!       call dm_view(A,ierr)
!       if(myrank==0) print *, ">B="
!       call dm_view(B,ierr)
!       if(myrank==0) print *, ">C=dm_setcol(A,0,B)"
!       call dm_view(C,ierr)
!       if(myrank==0) print *, ">D=dm_setcol(A,m-1,dm_m2n(10,10+m-1))"
!       call dm_view(D,ierr)
!       if(myrank==0) print *, ">E=dm_setcol(A,m-1,dm_m2n(10,10+m))"
!       call dm_view(E,ierr)
!       if(myrank==0) print *, ">F=dm_setcol(A,1,dm_zeros(m,1))"
!       call dm_view(F,ierr)
!       if(myrank==0) print *, ">G=dm_setcol(A,1,dm_ones(m,1))"
!       call dm_view(G,ierr)
!       if(myrank==0) print *, ">H=dm_setcol(H,1,dm_ones(m,1,.false.))"
!       if(myrank==0) call dm_view(H,ierr)
!	endif
! 	call dm_destroy(A,ierr)
! 	call dm_destroy(B,ierr)
! 	call dm_destroy(C,ierr)
! 	call dm_destroy(D,ierr)
! 	call dm_destroy(E,ierr)
! 	call dm_destroy(F,ierr)
! 	call dm_destroy(G,ierr)
! 	call dm_destroy(H,ierr)


!   if(myrank==0) print *, "==============Test dm_test================"
!   call dm_test(m,m,ierr)

	call dm_finalize(ierr)
end program
