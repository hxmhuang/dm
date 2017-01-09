
module dm_test
  use dm
  use dm_type
  use dm_op
  implicit none

#include "mat_type.h"
contains
  subroutine test_dm_zeros()
    type(Matrix)    		:: A,B,C,D
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)

    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_zeros==============="
    A=dm_zeros(m,n,k)
    B=dm_zeros(m,n,k,.true.)
    C=dm_zeros(m,n,2,.false.)
    D=dm_zeros(m,m,m)
    if(debug) then
       if(myrank==0) print *, ">A=dm_zeros(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_zeros(m,n,k,.true.)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_zeros(m,n,2,.false.)"
       if(myrank==0) call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_zeros(m,m,m)"
       call dm_view(D,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
  end subroutine test_dm_zeros


  subroutine test_dm_ones()
    type(Matrix)    		:: A,B,C,D,E
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_ones================"
    A=dm_ones(m,n,k)
    B=dm_ones(m,n,k,.true.)
    C=dm_ones(m,n,1,.false.)
    D=dm_ones(m,m,k)
    E=dm_ones(n,n,k)
    if(debug) then
       if(myrank==0) print *, ">A=dm_ones(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_ones(m,n,k,.true.)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_ones(m,n,2,.false.)"
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
  end subroutine test_dm_ones

  subroutine test_dm_eye()
    type(Matrix)    		:: A,B,C,D,E,F,G,H
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_eye================="
    A=dm_eye(m,n,k)	
    B=dm_eye(m,m*2,1)	
    C=dm_eye(2*m,m,1)	
    D=dm_eye(m,m*2,1,.true.)	
    E=dm_eye(2*m,m,1,.false.)	
    F=dm_eye(m,m*2,2)	
    G=dm_eye(2*m,m,3)
    H=dm_eye(m,m,1)	
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
       if(myrank==0) print *, ">H=dm_eye(m,m,1)"
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

  end subroutine test_dm_eye

  subroutine test_dm_copy()
    type(Matrix)    		:: A,B
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

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

  end subroutine test_dm_copy

  subroutine test_dm_add()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

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
    V=dm_eye(m,n,2,.false.)+dm_eye(m,n,2,.false.)
    W=2+dm_eye(m,n,2,.false.)
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
       if(myrank==0) print *, ">V=dm_eye(m,n,2,.false.)+dm_eye(m,n,2,.false.)"
       if(myrank==0) call dm_view(V,ierr)
       if(myrank==0) print *, ">W=2+dm_eye(m,n,2,.false.)"
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

  end subroutine test_dm_add

  subroutine test_dm_minus()
    type(Matrix)    		:: A,B,C,D,E,F,G,H
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

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
    V=dm_eye(m,n,2,.false.)-dm_eye(m,n,2,.false.)
    W=2-dm_eye(m,n,2,.false.)
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
       if(myrank==0) print *, ">V=dm_eye(m,n,2,.false.)-dm_eye(m,n,2,.false.)"
       if(myrank==0) call dm_view(V,ierr)
       if(myrank==0) print *, ">W=2-dm_eye(m,n,2,.false.)"
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

  end subroutine test_dm_minus

  subroutine test_dm_xjoin()
    type(Matrix)    		:: A,B,C,D,E,F,G,H
    type(Matrix)    		:: U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_xjoin==============="
    A=dm_eye(m,n,k)
    B=dm_eye(m+1,n,k)
    C=dm_eye(m-1,n,k)
    D=A .xj. B
    E=A .xj. C
    F=dm_eye(m,n,k) .xj. dm_eye(m,n,k)
    G=dm_eye(m,n,k) .xj. B
    H=A .xj. dm_eye(m,n,k)
    U=A .xj. A .xj. A
    V=A .xj. G
    W=dm_eye(m,n,2,.false.) .xj. dm_eye(m,n,2,.false.)
    if(debug) then
       if(myrank==0) print *, ">A=dm_eye(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_eye(m+1,n,k)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_eye(m-1,n,k)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=A .xj. B"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=A .xj. C"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=dm_eye(m,n,k) .xj. dm_eye(m,n,k)"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=dm_eye(m,n,k) .xj. B"
       call dm_view(G,ierr)
       if(myrank==0) print *, ">H=A .xj. dm_eye(m,n,k)"
       call dm_view(H,ierr)
       if(myrank==0) print *, ">U=A .xj. A .xj. A"
       call dm_view(U,ierr)
       if(myrank==0) print *, ">V=B .xj. G"
       call dm_view(V,ierr)
       if(myrank==0) print *, ">W=dm_eye(m,n,2,.false.) .xj. dm_eye(m,n,2,.false.)"
       if(myrank==0) call dm_view(W,ierr)
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
    call dm_destroy(V,ierr)
    call dm_destroy(W,ierr)

  end subroutine test_dm_xjoin

  subroutine test_dm_yjoin()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_yjoin==============="
    A=dm_eye(m,n,k)
    B=dm_eye(m,2*n,k)
    C=A .yj. B
    D=dm_eye(m,n,k) .yj. dm_eye(m,n,k)
    E=dm_eye(m,n,k) .yj. B
    F=A .yj. dm_eye(m,n,k)
    G=A .yj. A .yj. A
    H=B .yj. G
    U=dm_eye(m,n,2,.false.) .yj. dm_eye(m,n,2,.false.)
    if(debug) then
       if(myrank==0) print *, ">A=dm(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm(m,2*n,k)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=A .yj. B"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_eye(m,n,k) .yj. dm_eye(m,n,k)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=dm_eye(m,n,k) .yj. B"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=A .yj. dm_eye(m,n,k)"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=A .yj. A .yj. A"
       call dm_view(G,ierr)
       if(myrank==0) print *, ">H=B .yj. G"
       call dm_view(H,ierr)
       if(myrank==0) print *, ">U=dm_eye(m,n,2,.false.) .yj. dm_eye(m,n,2,.false.)"
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
  end subroutine test_dm_yjoin

  subroutine test_dm_zjoin()
    type(Matrix)    		:: A,B,C,D,E,F,G,H
    type(Matrix)    		:: U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_zjoin==============="
    A=dm_eye(m,n,k)
    B=dm_eye(m,n,k+1)
    C=dm_eye(m,n,k-1)
    D=A .zj. B
    E=A .zj. C
    F=dm_eye(m,n,k) .zj. dm_eye(m,n,k)
    G=dm_eye(m,n,k) .zj. B
    H=A .zj. dm_eye(m,n,k)
    U=A .zj. A .zj. A
    V=A .zj. G
    W=dm_eye(m,n,2,.false.) .zj. dm_eye(m,n,2,.false.)
    if(debug) then
       if(myrank==0) print *, ">A=dm_eye(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_eye(m,n,k+1)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_eye(m,n,k-1)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=A .zj. B"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=A .zj. C"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=dm_eye(m,n,k) .zj. dm_eye(m,n,k)"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=dm_eye(m,n,k) .zj. B"
       call dm_view(G,ierr)
       if(myrank==0) print *, ">H=A .zj. dm_eye(m,n,k)"
       call dm_view(H,ierr)
       if(myrank==0) print *, ">U=A .zj. A .zj. A"
       call dm_view(U,ierr)
       if(myrank==0) print *, ">V=B .zj. G"
       call dm_view(V,ierr)
       if(myrank==0) print *, ">W=dm_eye(m,n,2,.false.) .zj. dm_eye(m,n,2,.false.)"
       if(myrank==0) call dm_view(W,ierr)
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
    call dm_destroy(V,ierr)
    call dm_destroy(W,ierr)
  end subroutine test_dm_zjoin

  subroutine test_dm_mult()
    type(Matrix)    		:: A,B,C,D,E,F
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_mult================"
    A=dm_seqs(m,n,k)
    B=dm_seqs(n,n*2,k)
    C=A*B
    D=A*(dm_seqs(n,n*2,k))
    E=dm_seqs(m,n,k)*B
    F=dm_seqs(m,n,k)*dm_seqs(n,n*2,k) 
    X=A*2.0
    Y=2*A
    alpha=3.0
    Z=alpha*A
    U=A*alpha
    V=dm_seqs(m,n,2,.false.)*dm_seqs(n,n,2,.false.)
    W=3*dm_seqs(m,n,2,.false.)
    if(debug) then
       if(myrank==0) print *, ">A=dm_seqs(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_seqs(n,n*2,k)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=A*B"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=A*dm_seqs(n,n*2,k)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=dm_seqs(m,n,k)*B"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=dm_seqs(m,n,k)*dm_seqs(n,n*2,k)"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">X=A*2.0"
       call dm_view(X,ierr)
       if(myrank==0) print *, ">Y=2*A"
       call dm_view(Y,ierr)
       if(myrank==0) print *, ">Z=alpha*A"
       call dm_view(Y,ierr)
       if(myrank==0) print *, ">U=A*alpha"
       call dm_view(Y,ierr)
       if(myrank==0) print *, ">V=dm_seqs(m,n,2,.false.)*dm_seqs(n,n,2,.false.)"
       if(myrank==0) call dm_view(V,ierr)
       if(myrank==0) print *, ">W=3*dm_seqs(m,n,2,.false.)"
       if(myrank==0) call dm_view(W,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(X,ierr)
    call dm_destroy(Y,ierr)
    call dm_destroy(Z,ierr)
    call dm_destroy(U,ierr)
    call dm_destroy(V,ierr)
    call dm_destroy(W,ierr)


  end subroutine test_dm_mult

  subroutine test_dm_emult()
    type(Matrix)    		:: A,B,C,D,E,F,G,H
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_emult==============="
    A=dm_seqs(m,n,k)
    B=dm_seqs(m,n,k)*2
    C=A .em. B
    D=A .em. (dm_seqs(m,n,k))
    E=dm_seqs(m,n,k) .em. B
    F=dm_seqs(m,n,k) .em. dm_seqs(m,n,k) 
    G=A .em. A
    H=dm_seqs(m,n,2,.false.) .em. dm_seqs(m,n,2,.false.)
    if(debug) then
       if(myrank==0) print *, ">A=dm_seqs(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_seqs(m,n,k)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=A.*B"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=A.*dm_seqs(m,n,k)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=dm_seqs(m,n,k).*B"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=dm_seqs(m,n,k).*dm_seqs(m,n,k)"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=A.*A"
       call dm_view(G,ierr)
       if(myrank==0) print *, ">H=dm_seqs(m,n,2,.false.) .em. dm_seqs(m,n,2,.false.)"
       if(myrank==0) call dm_view(H,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_emult

  subroutine test_dm_ediv()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_ediv================"
    A=2*dm_eye(m,n,k)
    B=4*dm_eye(m,n,k)
    C=A .ed. B
    D=A .ed. (dm_eye(m,n,k))
    E=dm_eye(m,n,k) .ed. B
    F=dm_eye(m,n,k) .ed. dm_eye(m,n,k) 
    G=A .ed. A
    H=dm_eye(m,n,2,.false.) .ed. dm_eye(m,n,2,.false.)
    if(debug) then
       if(myrank==0) print *, ">A=2*dm_eye(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=4*dm_eye(m,n,k)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=A./B"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=A./dm_eye(m,n,k)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=dm_eye(m,n,k)./B"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=dm_eye(m,n,k)./dm_eye(m,n,k)"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=A./A"
       call dm_view(G,ierr)
       if(myrank==0) print *, ">H=dm_eye(m,n,2,.false.) .ed. dm_eye(m,n,2,.false.)"
       if(myrank==0) call dm_view(H,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_ediv

  subroutine test_dm_rep()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer                     :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_rep================="
    A=dm_seqs(m,n,k)
    B=dm_rep(A,1,2,1) 
    C=dm_rep(dm_seqs(m,n,k),2,3,1)
    D=dm_rep(dm_seqs(m,n,k),2,3,2)    
    E=dm_rep(dm_seqs(m,n,1,.false.),3,2,1) 
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_rep(A,1,2,1)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_rep(dm_seqs(m,n,k),2,3,1)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_rep(dm_seqs(m,n,k),2,3,2)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=dm_rep(dm_seqs(m,n,2,.false.),3,2,1)"
       if(myrank==0) call dm_view(E,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)    
  end subroutine test_dm_rep

  subroutine test_dm_seqs()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test dm_seqs================"
    A = dm_seqs(m, n, k)
    B = dm_seqs(m, n, k, .false.)
    if(debug) then
       if(myrank == 0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank == 0) print*, ">B="
       if(myrank == 0) call dm_view(B, ierr)
    end if
    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
  end subroutine test_dm_seqs

  subroutine test_dm_rand()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test dm_rand================"
    A = dm_rand(m, n, k)
    B = dm_rand(m, n, 2, .false.)
    if(debug) then
       if(myrank == 0) print*, ">A=dm_rand(m, n, k)"
       call dm_view(A, ierr)
       if(myrank == 0) print*, ">B=dm_rand(m, n, 2, .false.)"
       if(myrank == 0) call dm_view(B, ierr)
    end if
    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)

  end subroutine test_dm_rand

  subroutine test_dm_sum()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)
    real(kind=8) :: res
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test dm_sum================"
    A = dm_seqs(m, n, k)
    B = dm_seqs(m, n, 2, .false.)

    C = dm_sum(A, 1)
    D = dm_sum(A, 2)  
    res = dm_sum_all(A)
    
    if(debug) then
       if(myrank == 0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank == 0) print*, ">B="
       if(myrank == 0) call dm_view(B, ierr)
       if(myrank == 0) print*, ">C="
       call dm_view(C, ierr)
       if(myrank == 0) print*, ">D="
       call dm_view(D, ierr)
       if(myrank == 0) print*, ">dm_sum_all(A)=", res
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
    call dm_destroy(D, ierr)
  end subroutine test_dm_sum

  subroutine test_dm_axpy()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_axpy================"
    A=dm_seqs(m, n, k)
    B=dm_eye(m,n,k)
    C=dm_eye(m,n,k) 
    D=dm_eye(m,n,2,.false.) 

    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C="
       call dm_view(C,ierr)
       if(myrank==0) print*, ">D="
       if(myrank==0) call dm_view(D,ierr)
    endif

    call dm_axpy(B,2.0,A,ierr)
    if(debug) then
       if(myrank==0) print *, ">dm_axpy(B,2.0,A)"
       call dm_view(B,ierr)
    endif

    call dm_axpy(C,2,dm_seqs(m,n,k),ierr)
    if(debug) then
       if(myrank==0) print *, "dm_axpy(C,2,dm_seqs(m,n,k),ierr)"
       call dm_view(C,ierr)
    endif

    call dm_axpy(D,2,dm_seqs(m,n,2,.false.),ierr)
    if(debug) then
       if(myrank==0) print *, "dm_axpy(D,2,dm_seqs(m,n,2,.false.),ierr)"
       if(myrank==0) call dm_view(D,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
  end subroutine test_dm_axpy
  
  subroutine test_dm_aypx()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_aypx================"
    A=dm_seqs(m, n, k)
    B=dm_eye(m,n,k)
    C=dm_eye(m,n,k) 
    D=dm_eye(m,n,2,.false.) 

    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C="
       call dm_view(C,ierr)
       if(myrank==0) print*, ">D="
       if(myrank==0) call dm_view(D,ierr)
    endif

    call dm_aypx(B,2.0,A,ierr)
    if(debug) then
       if(myrank==0) print *, ">dm_aypx(B,2.0,A)"
       call dm_view(B,ierr)
    endif

    call dm_aypx(C,2,dm_seqs(m,n,k),ierr)
    if(debug) then
       if(myrank==0) print *, "dm_aypx(C,2,dm_seqs(m,n,k),ierr)"
       call dm_view(C,ierr)
    endif

    call dm_aypx(D,2,dm_seqs(m,n,2,.false.),ierr)
    if(debug) then
       if(myrank==0) print *, "dm_aypx(D,2,dm_seqs(m,n,2,.false.),ierr)"
       if(myrank==0) call dm_view(D,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
  end subroutine test_dm_aypx
  
  subroutine test_dm_trans()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_trans==============="
    A=dm_seqs(m,n,k)	
    B=dm_trans(A)
    C=dm_trans(dm_seqs(m,n,k))
    D=dm_trans(dm_seqs(m,n,2,.false.))
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B= dm_trans(A)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_trans(dm_seqs(m,n,k))"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_trans(dm_seqs(m,n,2,.false.))"
       if(myrank==0) call dm_view(D,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)

  end subroutine test_dm_trans
  
  subroutine test_dm_xyt()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_xyt================="
    A=dm_seqs(m,m,k)
    B=dm_ones(m,m,k)
    C=dm_xyt(A,B)
    D=dm_xyt(A,dm_ones(m,m,k))
    E=dm_xyt(B, dm_seqs(m,m,k))
    F=dm_xyt(dm_seqs(m,m,k),dm_ones(m,m,k))
    G=dm_xyt(A,A)
    H=dm_xyt(dm_seqs(m,m,2,.false.),dm_ones(m,m,2,.false.))

    if(debug) then
       if(myrank==0) print *, ">A=dm_seqs(m,m)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_ones(m,m)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_xyt(A,B)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_xyt(A,dm_ones(m,m))"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=dm_xyt(B, dm_seqs(m,m))"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=dm_xyt(dm_seqs(m,m),dm_ones(m,m))"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=dm_xyt(A,A)"
       call dm_view(G,ierr)
       if(myrank==0) print *, ">H=dm_xyt(dm_seqs(m,m,.false.),dm_ones(m,m,.false.))"
       if(myrank==0) call dm_view(H,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_xyt
  
  subroutine test_dm_xty()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_xty================="
    A=dm_seqs(m,m,k)
    B=dm_ones(m,m,k)
    C=dm_xty(A,B)
    D=dm_xty(A,dm_ones(m,m,k))
    E=dm_xty(dm_seqs(m,m,k),B)
    F=dm_xty(dm_seqs(m,m,k),dm_ones(m,m,k))
    G=dm_xty(A,A)
    H=dm_xty(dm_seqs(m,m,2,.false.),dm_ones(m,m,2,.false.))

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
       if(myrank==0) print *, ">H=dm_xty(dm_seqs(m,m,.false.),dm_ones(m,m,.false.))"
       if(myrank==0) call dm_view(H,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_xty


  subroutine test_dm_exp()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_exp================="
    A=dm_seqs(m,n,k)
    B=dm_exp(A)
    C=dm_exp(dm_seqs(m,n,k))
    D=dm_exp(dm_seqs(m,n,2,.false.))
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_exp(A)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_exp(dm_seqs(m,n,k))"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_exp(dm_seqs(m,n,2,.false.))"
       if(myrank==0) call dm_view(D,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)

  end subroutine test_dm_exp


  subroutine test_dm_log()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_log================="
    A=dm_seqs(m,m,k)
    B=dm_log(A)
    C=dm_log(dm_seqs(m,m,k))
    D=dm_log(dm_seqs(m,m,2,.false.))
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_log(A)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_log(dm_seqs(m,m,k))"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_log(dm_seqs(m,m,2,.false.))"
       if(myrank==0) call dm_view(D,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
  end subroutine test_dm_log

  subroutine test_dm_sqrt()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_sqrt================"
    A=dm_seqs(m,m,k)
    B=dm_sqrt(A)
    C=dm_sqrt(dm_seqs(m,m,k))
    D=dm_sqrt(dm_seqs(m,m,2,.false.))
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_sqrt(A)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_sqrt(dm_seqs(m,m,k))"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_sqrt(dm_seqs(m,m,2,.false.))"
       if(myrank==0) call dm_view(D,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
  end subroutine test_dm_sqrt
  
  subroutine test_dm_squ()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
  
    if(myrank==0) print *, "==============Test dm_squ================="
    A=dm_seqs(m,m,k)
    B=dm_squ(A)
    C=dm_squ(dm_seqs(m,m,k))
    D=dm_squ(dm_seqs(m,m,2,.false.))
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_squ(A)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_squ(dm_seqs(m,m,k))"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_squ(dm_seqs(m,m,2,.false.))"
       if(myrank==0) call dm_view(D,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
  end subroutine test_dm_squ

  subroutine test_dm_cube()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_cube================"
    A=dm_seqs(m,m,k)
    B=dm_cube(A)
    C=dm_cube(dm_seqs(m,m,k))
    D=dm_cube(dm_seqs(m,m,2,.false.))
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_cube(A)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_cube(dm_seqs(m,m,k))"
       if(myrank==0) print *, ">D=dm_cube(dm_seqs(m,m,2,.false.))"
       if(myrank==0) call dm_view(D,ierr)
       call dm_view(C,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
  end subroutine test_dm_cube
  
  subroutine test_dm_solve()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
  
    if(myrank==0) print *, "==============Test dm_solve==============="
    A=dm_rand(m,m,k, .true.)
    B=A * dm_ones(m,1,k, .true.)
    C=dm_solve(A,B)
    D=A .inv. B
    E=dm_rand(m,m,k, .true.)
    F=dm_rand(m,1,k, .true.)
    G=dm_solve(E,F)
    H=E .inv. F

    !call dm_view(dm_ones(m, 1, k, .true.), ierr)
    
    if(debug) then
       if(myrank==0) print *, ">A=dm_seqs(m,m,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_ones(m,1,k)"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_solve(A,B)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=A .inv. B"
       call dm_view(D,ierr)

       if(myrank==0) print *, ">E=dm_rand(m,m,2,.false.)"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=3*dm_rand(m,1,2,.false.)"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=dm_solve(E,F)"
       call dm_view(G,ierr)
       if(myrank==0) print *, ">H=E .inv. F"
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
  end subroutine test_dm_solve
  
  subroutine test_dm_setvalue()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
  
    if(myrank==0) print *, "==============Test dm_setvalue============"
    A=dm_eye(m,n,k)
    B=dm_eye(m,n,k)
    C=dm_eye(m,n,k)
    D=dm_eye(m,n,2,.false.)
    call dm_setvalue(A,1,1,0,5,ierr)
    call dm_setvalue(B,1,1,0,6.1,ierr)
    call dm_setvalue(C,1,1,0,real(7,kind=8),ierr)
    call dm_setvalue(D,1,0,1,8,ierr)
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C="
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D="
       if(myrank==0) call dm_view(D,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
  end subroutine test_dm_setvalue
  
  subroutine test_dm_getsub()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
  
    if(myrank==0) print *, "==============Test dm_getsub=============="
    A=dm_seqs(m,n,k)
    B=dm_getsub(A,(/0,1/),(/0,1/),(/0/))
    ! call dm_view(A, ierr)
    ! call dm_view(B, ierr)
    ! print*, B%nx, B%ny, B%nz
    C=dm_getsub(A,(/0,2,1/),(/1,0/),(/0,1/))
    D=dm_seqs(m,m,2,.false.)
    E=dm_getsub(D,(/0,1,2/),(/0,1/),(/0/))
    F=dm_getsub(D,(/0,2,1/),(/1,0/),(/0/))
    if(debug) then
       if(myrank==0) print *, ">A=dm_seqs(m,n,k)"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=dm_getsub(A,(/0,1/),(/0,1/),(/0/))"
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_getsub(A,(/0,2,1/),(/1,0/),(/0,1/))"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_seqs(m,m,2,.false.)"
       if(myrank==0) call dm_view(D,ierr)
       if(myrank==0) print *, ">E=dm_getsub(D,(/0,1,2/),(/0,1/),(/0/))"
       if(myrank==0) call dm_view(E,ierr)
       if(myrank==0) print *, ">F=dm_getsub(D,(/0,2,1/),(/1,0/),(/0/))"
       if(myrank==0) call dm_view(F,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)

  end subroutine test_dm_getsub

  subroutine test_dm_setvalues()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_setvalues==========="
    A=dm_ones(2*m,2*m, 2)
    B=dm_ones(2*m,2*m, 2, .false.)
    C=dm_ones(2*m,2*m, 2)

    call dm_setvalues(A,(/0,2/),(/1,2/),(/0/),(/9,8,7,6/),ierr)	
    call dm_setvalues(B,(/0,2/),(/1,2/),(/1/),(/9,8,7,6/),ierr)	
    call dm_setvalues(C,(/0,2/),(/1,2/),(/1,0/),(/9,8,7,6,0,1,2,3/),ierr)

    if(debug) then
       if(myrank==0) print *, &
            ">A=dm_setvalues(A,(/0,2/),(/1,2/),(/0/),(/9,8,7,6/),ierr)	"
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B=&
            &dm_setvalues(B,(/0,2/),(/1,2/),(/1/),(/9,8,7,6/),ierr)"
       if(myrank==0) call dm_view(B,ierr)
       if(myrank==0) print*, ">C=&
            &dm_setvalues(C,(/0,2/),(/1,2/),(/1,0/),(/9,8,7,6,0,1,2,3/),ierr)"
       call dm_view(C, ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)

  end subroutine test_dm_setvalues

  subroutine test_dm_getvalues()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_getvalues==========="
    A=dm_seqs(2*m,2*m,1)
    B=dm_seqs(2*m,2*m,2)
    C=dm_seqs(2*m,2*m,2, .false.)
    
    allocate(array(2))
    array=0

    call dm_getvalues(A,(/0/),(/1,2/),(/0/),array,ierr)	
    if(debug) then
       if(myrank==0) print *,">A=dm_seqs(2*m,2*m,1)"
       call dm_view(A,ierr)
       if(myrank==0) &
            print *, ">A=dm_getvalues(A,(/0/),(/1,2/),(/0/),array,ierr) : ",array
    endif

    deallocate(array)
    
    allocate(array(8))
    call dm_getvalues(B,(/0,1/),(/3,2/),(/0,1/),array,ierr)
    
    if(debug) then
       if(myrank==0) print *,">B=dm_seqs(2*m,2*m,2)" 
       call dm_view(B,ierr)
       if(myrank==0) &
            print *, ">B=dm_getvalues(B,(/0,1/),(/3,2/),(/0,1/),array,ierr) : ",array
    endif
    deallocate(array)
    
    allocate(array(8))
    call dm_getvalues(C,(/0,3/),(/3,2/),(/0,1/), array,ierr)
    
    if(debug) then
       if(myrank==0) print*, ">C=dm_seqs(2*m,2*m,3,.false.)"
       if(myrank==0) call dm_view(C,ierr)
       if(myrank==0) then
          print *, ">C=dm_getvalues(C,(/0,3/),(/3,2/),(/0,1/),array,ierr) : ", array
       endif
    endif
    deallocate(array)
    
    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)

  end subroutine test_dm_getvalues

  subroutine test_dm_norm()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_norm================"
    A=dm_seqs(m, m, 2)
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
    B=dm_seqs(m,m,1,.false.)
    a1=dm_norm_1(B)
    a2=dm_norm_2(B)
    a3=dm_norm_inf(B) 
    if(debug) then
       if(myrank==0) print *, ">B="
       if(myrank==0) call dm_view(B,ierr)
       if(myrank==0) print *, ">NORM_1=",a1
       if(myrank==0) print *, ">NORM_2=",a2
       if(myrank==0) print *, ">NORM_INF=",a3
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
  end subroutine test_dm_norm


  subroutine test_dm_lt()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_lt=================="
    A=dm_seqs(m,m,2)
    B=5*dm_ones(m,m,2)
    C=(A<B)
    D=(A<5)
    E=(A<5.0)
    F=(A<real(5.0, kind=8))
    G=(dm_seqs(m, m, 1, .false.) < 3.0*dm_ones(m, m, 1, .false.))
    H=(dm_seqs(m, m, 2, .false.) < 5.0)

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
       if(myrank==0) print *, &
            ">G=(dm_seqs(m,m,1.false.) < 3.0*dm_ones(m,m,1,.false.))"
       if(myrank==0) call dm_view(G,ierr)
       if(myrank==0) print *, ">H=(dm_seqs(m,m,2,.false.) < 5)"
       if(myrank==0) call dm_view(H,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_lt

  subroutine test_dm_le()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_le=================="
    A=dm_seqs(m,n,1)
    B=5*dm_ones(m,n,1)
    C=(A<=B)
    D=(A<=4)
    E=(A<=4.0)
    F=(A<=real(5.0,kind=8))
    G=(dm_seqs(m,n,1,.false.) <= 5.0*dm_seqs(m,n,1,.false.))
    H=(dm_seqs(m,n,2,.false.) <= 5)

    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=(A<=B)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=(A<=4)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=(A<=4.0)"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=(A<=real(5.0,kind=8))"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=(dm_seqs(m,n,1.false.) <= 5.0*dm_seqs(m,n,1,.false.))"
       if(myrank==0) call dm_view(G,ierr)
       if(myrank==0) print *, ">H=(dm_seqs(m,n,2.false.) <= 5)"
       if(myrank==0) call dm_view(H,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_le

  subroutine test_dm_gt()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank==0) print *, "==============Test dm_gt=================="

    A=dm_seqs(m,n,1)
    B=4*dm_ones(m,n,1)
    C=(A>B)
    D=(A>3)
    E=(A>4.0)
    F=(A>real(4.0,kind=8))
    G=(dm_seqs(m,n,1,.false.) > 4.0*dm_ones(m,n,1,.false.))
    H=(dm_seqs(m,n,2,.false.) > 4)

    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=(A>B)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=(A>3)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=(A>4.0)"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=(A>real(4.0,kind=8))"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=(dm_seqs(m,n,1,.false.) > 4.0*(m,n,1,.false.))"
       if(myrank==0) call dm_view(G,ierr)
       if(myrank==0) print *, ">H=(dm_seqs(m,n,2,.false.) > 4)"
       if(myrank==0) call dm_view(H,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_gt

  subroutine test_dm_ge()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_ge=================="
    A=dm_seqs(m,n,1)
    B=5*dm_ones(m,n,1)
    C=(A>=B)
    D=(A>=3)
    E=(A>=4.0)
    F=(A>=real(4.0,kind=8))
    G=(dm_seqs(m,n,1,.false.) >= 4.0*dm_ones(m,n,1,.false.))
    H=(dm_seqs(m,n,1,.false.) >= 5)

    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=(A>=B)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=(A>=3)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=(A>=4.0)"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=(A>=real(4.0,kind=8))"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=(dm_seqs(m,n,1,.false.) >= 4.0*dm_ones(m,n,1,.false.))"
       if(myrank==0) call dm_view(G,ierr)
       if(myrank==0) print *, ">H=(dm_seqs(m,n,1,.false.) >= 5)"
       if(myrank==0) call dm_view(H,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_ge

  subroutine test_dm_eq()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_eq=================="
    A=dm_seqs(m,n,1)
    B=5*dm_ones(m,n,1)
    C=(A==B)
    D=(A==3)
    E=(A==3.0)
    F=(A==real(3.0,kind=8))
    G=(dm_seqs(m,n,1,.false.) == 4.0*dm_ones(m,n,1,.false.))
    H=(dm_seqs(m,n,2,.false.) == 7)
    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=(A==B)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=(A==3)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=(A==3.0)"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=(A==real(3.0,kind=8))"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=(dm_seqs(m,n,1,.false.)==4.0*dm_ones(m,n,1,.false.))"
       if(myrank==0) call dm_view(G,ierr)
       if(myrank==0) print *, ">H=(dm_seqs(m,n,2,.false.) == 7)"
       if(myrank==0) call dm_view(H,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
  end subroutine test_dm_eq

  subroutine test_dm_nq()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_nq=================="
    A=dm_seqs(m,n,1)
    B=5*dm_ones(m,n,1)
    C=(A/=B)
    D=(A/=4)
    E=(A/=5.0)
    F=(A/=real(5.0,kind=8))
    G=(dm_seqs(m,n,1,.false.) /= 2.0*dm_ones(m,n,1,.false.))
    H=(dm_seqs(m,n,2,.false.) /= 7)

    if(debug) then
       if(myrank==0) print *, ">A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=(A/=B)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=(A/=4)"
       call dm_view(D,ierr)
       if(myrank==0) print *, ">E=(A/=5.0)"
       call dm_view(E,ierr)
       if(myrank==0) print *, ">F=(A/=real(5.0,kind=8))"
       call dm_view(F,ierr)
       if(myrank==0) print *, ">G=(dm_seqs(m,n,1,.false.)/=2.0*dm_ones(m,n,1,.false.))"
       if(myrank==0) call dm_view(G,ierr)
       if(myrank==0) print *, ">H=(dm_seqs(m,n,1,.false.)/=7)"
       if(myrank==0) call dm_view(H,ierr)
    endif
    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)

  end subroutine test_dm_nq

  subroutine test_dm_max_min()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer                     :: pos1(3), pos2(3), pos3(3), pos4(3)
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)
    real(kind=8) :: m1, m2, m3, m4
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_max_min=============="

    A = dm_rand(m, n, k, .true.)
    B = dm_seqs(m, n, k, .false.)
    
    call dm_max(A, m1, pos1, ierr)
    call dm_min(A, m2, pos2, ierr)
    call dm_max(B, m3, pos3, ierr)
    call dm_min(B, m4, pos4, ierr)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)

       if(myrank==0) then
          print*, ">max(A)=", m1, " POS=", pos1
          print*, ">min(A)=", m2, " POS=", pos2          
       endif

       if(myrank==0) then
          print*, ">B="          
          call dm_view(B, ierr)
          print*, ">max(B)=", m3, " POS=", pos3
          print*, ">min(B)=", m4, " POS=", pos4          
       endif
    endif
    
    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
  end subroutine test_dm_max_min

  subroutine test_dm_sparse()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank==0) print *, "==============Test dm_sparse=============="

    A=dm_seqs(2,1,1) .xj. (dm_seqs(3,1,1)+2)
    B=(dm_seqs(2,1,1)) .xj. (dm_seqs(3,1,1)+3)
    C=(dm_zeros(5,1,1))
    D=(dm_seqs(5,1,1))
    E=dm_sparse(A,B,C,D,6,6,1)

    F=dm_seqs(3,1,1,.false.) .xj. (dm_seqs(3,1,1,.false.)+2)
    G=dm_seqs(3,1,1,.false.) .xj. (dm_seqs(3,1,1,.false.)+4)
    H=dm_zeros(6,1,1,.false.) + 1
    II=dm_seqs(6,1,1,.false.)
    KK=dm_sparse(F,G,H,II,8,8,2)

    if(debug) then
       if(myrank==0) print *, ">Ind_i:A="
       call dm_view(A,ierr)
       if(myrank==0) print *, ">Ind_j:B="
       call dm_view(B,ierr)
       if(myrank==0) print *, ">C=dm_zeros(1,5,1)"
       call dm_view(C,ierr)
       if(myrank==0) print *, ">D=dm_seqs(1,5,1)"
       call dm_view(D, ierr)
       if(myrank==0) print *, ">E=dm_sparse(1,5,1)"
       call dm_view(E, ierr)

       if(myrank==0) print *, ">Ind_i:F="
       if(myrank==0) call dm_view(F,ierr)
       if(myrank==0) print *, ">Ind_j:G="
       if(myrank==0) call dm_view(G,ierr)
       if(myrank==0) print *, ">Ind_k:H"
       if(myrank==0) call dm_view(H,ierr)
       if(myrank==0) print *, ">Ind_A:I"
       if(myrank==0) call dm_view(II,ierr)
       if(myrank==0) print *, ">Ind_B:K"
       if(myrank==0) call dm_view(KK,ierr)
    endif

    call dm_destroy(A,ierr)
    call dm_destroy(B,ierr)
    call dm_destroy(C,ierr)
    call dm_destroy(D,ierr)
    call dm_destroy(E,ierr)
    call dm_destroy(F,ierr)
    call dm_destroy(G,ierr)
    call dm_destroy(H,ierr)
    call dm_destroy(II,ierr)
    call dm_destroy(KK,ierr)
  end subroutine test_dm_sparse

  subroutine test_dm_save()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)
    character(len=12) :: file_name
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank == 0) print*, "==============Test dm_save================"
    A = dm_seqs(4,3,3)
    B = dm_ones(m, n, k)

    if(myrank==0 .and. debug) print*, ">>writing to file Test_A.nc"  
    call dm_save("Test_A.nc", "wind", A, ierr)

    if(myrank==0 .and. debug) print*, ">>writing to file Test_B.nc"    
    call dm_save("Test_B.nc", "wind", B, ierr)

    write(file_name, "(A,I3.3,A)") "Test_C", myrank, ".nc"
    
    if(myrank==0 .and. debug) print*, ">>writing to file ", file_name
    call dm_save(file_name, "wind", dm_seqs(m, n, k, .false.), ierr)

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
  end subroutine test_dm_save

  subroutine test_dm_load()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)
    character(len=12) :: file_name
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test dm_load================"

    call dm_load("Test_A.nc", "wind", A, .true.,  ierr)
    call dm_load("Test_B.nc", "wind", B, .true.,  ierr)

    write(file_name, "(A,I3.3,A)") "Test_C", myrank, ".nc"
    call dm_load(file_name, "wind", C, .false., ierr)

    if(debug) then
       if(myrank==0) print '(A, I2, I2, I2)', ">dim(A)=", A%nx, A%ny, A%nz
       call dm_view(A, ierr)

       if(myrank==0) print '(A, I2, I2, I2)', ">dim(B)=", B%nx, B%ny, B%nz
       call dm_view(B, ierr)

       if(myrank==0) print '(A, I2, I2, I2)', ">dim(C)=", C%nx, C%ny, C%nz     
       if(myrank==0) call dm_view(C, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
  end subroutine test_dm_load

  subroutine test_dm_save3d()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)
    character(len=12) :: file_name
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
   
    if(myrank == 0) print*, "==============Test dm_save3d=============="

    E = dm_seqs(m, n, k)
    F = dm_seqs(m*2,n+1,k)
    G = dm_seqs(m+1,n*2,k+1)
    H = dm_seqs(m,n,k*2)

    if(myrank==0 .and. debug) print*, ">>writing data into Test_E.nc"
    call dm_save3d("Test_E.nc", "wind", E, ierr)

    if(myrank==0 .and. debug) print*, ">>writing data into Test_F.nc"
    call dm_save3d("Test_F.nc", "wind", F, ierr)

    if(myrank==0 .and. debug) print*, ">>writing data into Test_G.nc"
    call dm_save3d("Test_G.nc", "wind", G, ierr)

    if(myrank==0 .and. debug) print*, ">>writing data into Test_H.nc"  
    call dm_save3d("Test_H.nc", "wind", H, ierr)

    write(file_name, "(A,I3.3,A)") "Test_I", myrank, ".nc"
    if(debug) print*, ">>writing data into ",file_name
    call dm_save3d(file_name, "wind", dm_seqs(m,n,k, .false.), ierr)

    call dm_destroy(E, ierr)
    call dm_destroy(F, ierr)
    call dm_destroy(G, ierr)
    call dm_destroy(H, ierr)
  end subroutine test_dm_save3d

  subroutine test_dm_load3d()
    type(Matrix)    		:: A,B,C,D,E,F,G,H,II
    type(Matrix)    		:: X,Y,Z,U,V,W 
    integer         		:: m,n,k
    integer                     :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)
    character(len=12) :: file_name
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test dm_load3d=============="

    call dm_load3d("Test_E.nc", "wind", E, .true., ierr)
    call dm_load3d("Test_F.nc", "wind", F, .true., ierr)
    call dm_load3d("Test_G.nc", "wind", G, .true., ierr)
    write(file_name, "(A,I3.3,A)") "Test_I",myrank,".nc"
    call dm_load3d(file_name, "wind", II, .false., ierr)

    if(debug) then
       if(myrank==0) print "(A,I3,A,I3,A,I3)", ">dim(E) =", E%nx,",",E%ny,",",E%nz     
       call dm_view(E, ierr)
       if(myrank==0) print "(A,I3,A,I3,A,I3)", ">dim(F) =", F%nx,",",F%ny,",",F%nz
       call dm_view(F, ierr)
       if(myrank==0) print "(A,I3,A,I3,A,I3)", ">dim(G) =", G%nx,",",G%ny,",",G%nz
       call dm_view(G, ierr)
       if(myrank==0) print "(A,I3,A,I3,A,I3)", ">dim(II) =", II%nx,",",II%ny,",",II%nz
       if(myrank==0) call dm_view(II, ierr)     
    endif

    call dm_destroy(E, ierr)
    call dm_destroy(F, ierr)
    call dm_destroy(G, ierr)
    call dm_destroy(II, ierr)
  end subroutine test_dm_load3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The following test the operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine test_OP_AXF()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_AXF=============="

    A = OP_AXF(m*2+1, m*2+1, k)
    call dm_save3d("operators.nc", "AXF", A, ierr)
    call dm_load3d("operators.nc", "AXF", B,.true., ierr)
    C = OP_AXF(m*2+1, m*2+1, k, .false.)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="       
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       if(myrank==0) call dm_view(C, ierr)
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
  end subroutine 

  
  subroutine test_OP_AYF()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_AYF=============="

    A = OP_AYF(m*2+1, n*2+1, k)
    call dm_save3d("operators.nc", "AYF", A, ierr)
    call dm_load3d("operators.nc", "AYF", B, .true., ierr)
    C = OP_AYF(m*2+1, n*2, k, .false.)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="       
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       if(myrank==0) call dm_view(C, ierr)
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
  end subroutine 

  subroutine test_OP_AZF()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_AZF=============="

    ! A = OP_AXF(m*2, m*2, k)
    ! call dm_save3d("operators.nc", "AXF", A, ierr)
    ! call dm_load3d("operators.nc", "AXF", B,.true., ierr)
    ! C = OP_AXF(m*2, m*2, k, .false.)
    
    ! if(debug) then
    !    if(myrank==0) print*, ">A="
    !    call dm_view(A, ierr)
    !    if(myrank==0) print*, ">B="       
    !    call dm_view(B, ierr)
    !    if(myrank==0) print*, ">C="
    !    if(myrank==0) call dm_view(C, ierr)
    ! end if

    ! call dm_destroy(A, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(C, ierr)
    
  end subroutine 

  subroutine test_OP_AXB()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_AXB=============="

    A = OP_AXB(m*2+1, n*2+1, k)
    call dm_save3d("operators.nc", "AXB", A, ierr)
    call dm_load3d("operators.nc", "AXB", B, .true., ierr)
    C = OP_AXB(m*2+1, n*2+1, k, .false.)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="       
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       if(myrank==0) call dm_view(C, ierr)
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
  end subroutine 

  subroutine test_OP_AYB()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_AYB=============="

    A = OP_AYB(m*2+1, m*2+1, k)
    call dm_save3d("operators.nc", "AYB", A, ierr)
    call dm_load3d("operators.nc", "AYB", B, .true., ierr)
    C = OP_AYB(m*2+1, n*2+1, k, .false.)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="       
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       if(myrank==0) call dm_view(C, ierr)
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
  end subroutine 

  subroutine test_OP_AZB()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_AZB=============="

    ! A = OP_AXF(m*2, m*2, k)
    ! call dm_save3d("operators.nc", "AXF", A, ierr)
    ! call dm_load3d("operators.nc", "AXF", B,.true., ierr)
    ! C = OP_AXF(m*2, m*2, k, .false.)
    
    ! if(debug) then
    !    if(myrank==0) print*, ">A="
    !    call dm_view(A, ierr)
    !    if(myrank==0) print*, ">B="       
    !    call dm_view(B, ierr)
    !    if(myrank==0) print*, ">C="
    !    if(myrank==0) call dm_view(C, ierr)
    ! end if

    ! call dm_destroy(A, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(C, ierr)
    
  end subroutine 


  subroutine test_OP_DXB()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_DXB=============="

    A = OP_DXB(m*2+1, n*2+1, k)
    call dm_save3d("operators.nc", "DXB", A, ierr)
    call dm_load3d("operators.nc", "DXB", B,.true., ierr)
    C = OP_DXB(m*2+1, n*2+1, k, .false.)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="       
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       if(myrank==0) call dm_view(C, ierr)
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
  end subroutine 

  subroutine test_OP_DYB()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_DYB=============="

    A = OP_DYB(m*2+1, m*2+1, k)
    call dm_save3d("operators.nc", "DYB", A, ierr)
    call dm_load3d("operators.nc", "DYB", B,.true., ierr)
    C = OP_DYB(m*2+1, n*2+1, k, .false.)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="       
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       if(myrank==0) call dm_view(C, ierr)
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
  end subroutine 


  subroutine test_OP_DZB()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_AZB=============="

    ! A = OP_AXF(m*2, m*2, k)
    ! call dm_save3d("operators.nc", "AXF", A, ierr)
    ! call dm_load3d("operators.nc", "AXF", B,.true., ierr)
    ! C = OP_AXF(m*2, m*2, k, .false.)
    
    ! if(debug) then
    !    if(myrank==0) print*, ">A="
    !    call dm_view(A, ierr)
    !    if(myrank==0) print*, ">B="       
    !    call dm_view(B, ierr)
    !    if(myrank==0) print*, ">C="
    !    if(myrank==0) call dm_view(C, ierr)
    ! end if

    ! call dm_destroy(A, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(C, ierr)
    
  end subroutine 

  subroutine test_OP_DXF()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_DXF=============="

    A = OP_DXF(m*2+1, n*2+1, k)
    call dm_save3d("operators.nc", "DXF", A, ierr)
    call dm_load3d("operators.nc", "DXF", B, .true., ierr)
    C = OP_DXF(m*2+1, n*2+1, k, .false.)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="       
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       if(myrank==0) call dm_view(C, ierr)
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
  end subroutine 

  subroutine test_OP_DYF()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_DYF=============="

    A = OP_DYF(m*2+1, n*2+1, k)
    call dm_save3d("operators.nc", "DYF", A, ierr)
    call dm_load3d("operators.nc", "DYF", B,.true., ierr)
    C = OP_DYF(m*2+1, n*2+1, k, .false.)
    
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="       
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       if(myrank==0) call dm_view(C, ierr)
    end if

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
    
  end subroutine 

  subroutine test_OP_DZF()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_DZF=============="

    ! A = OP_AXF(m*2, m*2, k)
    ! call dm_save3d("operators.nc", "AXF", A, ierr)
    ! call dm_load3d("operators.nc", "AXF", B,.true., ierr)
    ! C = OP_AXF(m*2, m*2, k, .false.)
    
    ! if(debug) then
    !    if(myrank==0) print*, ">A="
    !    call dm_view(A, ierr)
    !    if(myrank==0) print*, ">B="       
    !    call dm_view(B, ierr)
    !    if(myrank==0) print*, ">C="
    !    if(myrank==0) call dm_view(C, ierr)
    ! end if

    ! call dm_destroy(A, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(C, ierr)
    
  end subroutine 


  subroutine test_OP_DXC()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_DXC=============="

    ! A = OP_AXF(m*2, m*2, k)
    ! call dm_save3d("operators.nc", "AXF", A, ierr)
    ! call dm_load3d("operators.nc", "AXF", B,.true., ierr)
    ! C = OP_AXF(m*2, m*2, k, .false.)
    
    ! if(debug) then
    !    if(myrank==0) print*, ">A="
    !    call dm_view(A, ierr)
    !    if(myrank==0) print*, ">B="       
    !    call dm_view(B, ierr)
    !    if(myrank==0) print*, ">C="
    !    if(myrank==0) call dm_view(C, ierr)
    ! end if

    ! call dm_destroy(A, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(C, ierr)
    
  end subroutine 

  subroutine test_OP_DYC()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_DYC=============="

    ! A = OP_AXF(m*2, m*2, k)
    ! call dm_save3d("operators.nc", "AXF", A, ierr)
    ! call dm_load3d("operators.nc", "AXF", B,.true., ierr)
    ! C = OP_AXF(m*2, m*2, k, .false.)
    
    ! if(debug) then
    !    if(myrank==0) print*, ">A="
    !    call dm_view(A, ierr)
    !    if(myrank==0) print*, ">B="       
    !    call dm_view(B, ierr)
    !    if(myrank==0) print*, ">C="
    !    if(myrank==0) call dm_view(C, ierr)
    ! end if

    ! call dm_destroy(A, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(C, ierr)
    
  end subroutine 

  subroutine test_OP_DZC()
    type(Matrix)    		:: A,B,C
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    real(kind=8)    		:: a1,a2,a3
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable  :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test OP_DZC=============="

    ! A = OP_AXF(m*2, m*2, k)
    ! call dm_save3d("operators.nc", "AXF", A, ierr)
    ! call dm_load3d("operators.nc", "AXF", B,.true., ierr)
    ! C = OP_AXF(m*2, m*2, k, .false.)
    
    ! if(debug) then
    !    if(myrank==0) print*, ">A="
    !    call dm_view(A, ierr)
    !    if(myrank==0) print*, ">B="       
    !    call dm_view(B, ierr)
    !    if(myrank==0) print*, ">C="
    !    if(myrank==0) call dm_view(C, ierr)
    ! end if

    ! call dm_destroy(A, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(C, ierr)
    
  end subroutine 

  subroutine test_AXF()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test AXF=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = AXF(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AXF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_AXB()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test AXB=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = AXB(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AXF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 
  
  subroutine test_AYF()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test AYF=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = AYF(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_AYB()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test AYB=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = AYB(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_AZF()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test AZF=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = AZF(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_AZB()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test AZB=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = AZB(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  
  subroutine test_DXF()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DXF=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DXF(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  
  subroutine test_DXB()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DXB=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DXB(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_DXC()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DXC=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DXC(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 
  
  subroutine test_DYF()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DYF=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DYF(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_DYB()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DYB=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DYB(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_DYC()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DYC=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DYC(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_DZF()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DZF=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DZF(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_DZB()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)

    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DZB=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DZB(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_DZC()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1,B1,C1    
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)
    integer :: i
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test DZC=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = DZC(A)

    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_CSUM()
    use dm_op
    type(Matrix)    		:: A, B, C
    type(Matrix)    		:: A1, A2, A3
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)
    integer :: i
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test CSUM=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    A1 = CSUM(A, 1) !type 1, corresponding to SUM1 in matlab code
    A2 = CSUM(A, 2) !type 2, corresponding to SUM2 in matlab code
    A3 = CSUM(A, 3) !type 3
    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">A1="
       call dm_view(A1, ierr)
       if(myrank==0) print*, ">A2="
       call dm_view(A2, ierr)
       if(myrank==0) print*, ">A3="
       call dm_view(A3, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(A1, ierr)
    call dm_destroy(A2, ierr)
    call dm_destroy(A3, ierr)
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 

  subroutine test_SHIFT()
    use dm_op
    type(Matrix)    		:: A, B, C, D, E, F, G
    type(Matrix)    		:: A1, A2, A3
    integer         		:: m,n,k
    integer :: myrank, mysize
    real(kind=8)    		:: ep,alpha
    logical         		:: debug = .false.
    integer         		:: ierr
    real(kind=8),allocatable    :: array(:)
    integer :: i
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-k',k,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)

    if(myrank == 0) print*, "==============Test SHIFT=============="

    A = dm_seqs(2*m+1, 2*n+1, k)
    B = SHIFT(A, 1,  1) !shift on x-axis, forward
    C = SHIFT(A, 1, -1) !shift on x-axis, backward
    D = SHIFT(A, 2,  1) !shift on y-axis, forward
    E = SHIFT(A, 2, -1) !shift on y-axis, backward
    F = SHIFT(A, 3,  1) !shift on z-axis, forward
    G = SHIFT(A, 3, -1) !shift on z-axis, backward
    
    ! B = dm_seqs(2*m+1, 2*n+1, k, .false.)
    ! B1 = AYF(B)
    if(debug) then
       if(myrank==0) print*, ">A="
       call dm_view(A, ierr)
       if(myrank==0) print*, ">B="
       call dm_view(B, ierr)
       if(myrank==0) print*, ">C="
       call dm_view(C, ierr)
       if(myrank==0) print*, ">D="
       call dm_view(D, ierr)
       if(myrank==0) print*, ">E="
       call dm_view(E, ierr)
       if(myrank==0) print*, ">F="
       call dm_view(F, ierr)
       if(myrank==0) print*, ">G="
       call dm_view(G, ierr)
       ! if(myrank==0) print*, ">B="
       ! if(myrank==0) call dm_view(B, ierr)
       ! if(myrank==0) print*, ">B1="
       ! if(myrank==0) call dm_view(B1, ierr)
    endif

    call dm_destroy(A, ierr)
    call dm_destroy(B, ierr)
    call dm_destroy(C, ierr)
    call dm_destroy(D, ierr)
    call dm_destroy(E, ierr)
    call dm_destroy(F, ierr)
    call dm_destroy(G, ierr)        
    ! call dm_destroy(B, ierr)
    ! call dm_destroy(B1, ierr)
  end subroutine 
  
end module dm_test
