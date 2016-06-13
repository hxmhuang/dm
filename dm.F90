! -----------------------------------------------------------------------
! Distributed Matrix Wrapper 
! -----------------------------------------------------------------------
module dm 
    use dm_type
    use dm_mat
    implicit none
#include "mat_type.h"

    interface operator (+)
        module procedure dm_add1
        module procedure dm_add2
        module procedure dm_add3
        module procedure dm_add4
        module procedure dm_add5
        module procedure dm_add6
        module procedure dm_add7
    end interface

    interface operator (-)
        module procedure dm_minus1
        module procedure dm_minus2
        module procedure dm_minus3
        module procedure dm_minus4
        module procedure dm_minus5
        module procedure dm_minus6
        module procedure dm_minus7
    end interface

    interface operator (*)
        module procedure dm_mult1
        module procedure dm_mult2
        module procedure dm_mult3
        module procedure dm_mult4
        module procedure dm_mult5
        module procedure dm_mult6
        module procedure dm_mult7
    end interface
 
    interface operator (<)
        module procedure dm_lt1
        module procedure dm_lt2
        module procedure dm_lt3
        module procedure dm_lt4
    end interface
  
    interface operator (<=)
        module procedure dm_le1
        module procedure dm_le2
        module procedure dm_le3
        module procedure dm_le4
    end interface
    
    interface operator (>)
        module procedure dm_gt1
        module procedure dm_gt2
        module procedure dm_gt3
        module procedure dm_gt4
    end interface

    interface operator (>=)
        module procedure dm_ge1
        module procedure dm_ge2
        module procedure dm_ge3
        module procedure dm_ge4
    end interface
 
    interface operator (==)
        module procedure dm_eq1
        module procedure dm_eq2
        module procedure dm_eq3
        module procedure dm_eq4
    end interface
  
    interface operator (/=)
        module procedure dm_nq1
        module procedure dm_nq2
        module procedure dm_nq3
        module procedure dm_eq4
    end interface
  
	! element multiple
    interface operator (.to.)
        module procedure dm_to
    end interface
 
  	! element multiple
    interface operator (.em.)
        module procedure dm_emult
    end interface
     
    ! element divide 
    interface operator (.ed.)
        module procedure dm_ediv
    end interface
    
    ! join horizontally
    interface operator (.hj.)
        module procedure dm_hjoin
    end interface
    
	! join vertically
    interface operator (.vj.)
        module procedure dm_vjoin
    end interface
	! INV(A)*B
    interface operator (.inv.)
        module procedure dm_solve
    end interface

    interface dm_axpy
        module procedure dm_axpy1
        module procedure dm_axpy2
        module procedure dm_axpy3
    end interface

    interface dm_aypx
        module procedure dm_aypx1
        module procedure dm_aypx2
        module procedure dm_aypx3
    end interface

	interface dm_setvalue
        module procedure dm_setvalue1
        module procedure dm_setvalue2
        module procedure dm_setvalue3
    end interface

	interface dm_setvalues
        module procedure dm_setvalues1
        module procedure dm_setvalues2
        module procedure dm_setvalues3
    end interface
	
    interface assignment(=)
        module procedure dm_copy
    end interface
	
	interface operator (.to.)
        module procedure dm_add1 
    end interface



contains

! -----------------------------------------------------------------------
! Initialize the distributed matrix environment 
! -----------------------------------------------------------------------
subroutine dm_init(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)  ::  ierr 
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
end subroutine 


! -----------------------------------------------------------------------
! Get the rank number of the current process in the commmunicator 
! -----------------------------------------------------------------------
subroutine dm_comm_rank(myrank,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)         ::  myrank
    integer,intent(out)  		::  ierr 
	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
end subroutine 


! -----------------------------------------------------------------------
! Get the size of processes in the commmunicator 
! -----------------------------------------------------------------------
subroutine dm_comm_size(mysize,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)         ::  mysize
    integer,intent(out)  		::  ierr 
	call MPI_Comm_size(PETSC_COMM_WORLD,mysize,ierr)
end subroutine 


! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
subroutine dm_option_int(str,input,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*),intent(in) ::  str
    integer,intent(out)			::  input 
    integer,intent(out)			::  ierr 
    input=0
	call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,input,PETSC_NULL_BOOL,ierr)
end subroutine 


! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
subroutine dm_option_bool(str,input,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*),intent(in) ::  str
    logical,intent(out)			::  input 
    integer,intent(out)			::  ierr 
    input=.false.
	call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-debug',input,PETSC_NULL_BOOL,ierr)
end subroutine 


! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
subroutine dm_option_real(str,input,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*),intent(in) ::  str
    real(kind=8),intent(out)	::  input 
    integer,intent(out)			::  ierr 
    input=0.0
	call PetscOptionsGetReal(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,input,PETSC_NULL_BOOL,ierr)
end subroutine 


! -----------------------------------------------------------------------
! Finalize the distributed matrix environment 
! -----------------------------------------------------------------------
subroutine dm_finalize(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)			::  ierr 
    call PetscFinalize(ierr)
end subroutine 


! -----------------------------------------------------------------------
! Destroy a matrix to free the memory
! -----------------------------------------------------------------------
subroutine dm_destroy(A,ierr)
    implicit none
    type(Matrix),   intent(in)  ::  A
    integer,intent(out)			::  ierr 
    call mat_destroy(A%x,ierr)
end subroutine 


! -----------------------------------------------------------------------
! Print a matrix on screen
! -----------------------------------------------------------------------
subroutine dm_view(A,ierr) 
    implicit none
    type(Matrix),  intent(in)   ::  A
    integer,	   intent(out)	::  ierr 
    call mat_view(A%x,ierr)
end subroutine 

! -----------------------------------------------------------------------
! A=0 
! -----------------------------------------------------------------------
function dm_zeros(m,n) result(A)
    implicit none
    integer, 	intent(in)  ::  m,n 
    type(Matrix)            ::  A
    integer					::  ierr
    
    call mat_zeros(A%x,m,n,ierr)
    call dm_set_implicit(A,ierr)
end function


function dm_ones(m,n) result(A)
    implicit none
    integer,    intent(in)  ::  m,n
    type(Matrix)            ::  A
    integer					::  ierr
    
    call mat_ones(A%x,m,n,ierr)
    call dm_set_implicit(A,ierr)
end function


! -----------------------------------------------------------------------
! A=[0 1 2], This function is only used to generate the test data.
!   [3 4 5]
!   [6 7 8]
! -----------------------------------------------------------------------
function dm_seqs(m,n) result(A)
    implicit none
    integer,   intent(in)  	::  m,n
    type(Matrix)           	::  A
    integer					::  ierr

    call mat_seqs(A%x,m,n,ierr)
    call dm_set_implicit(A,ierr)
end function


! -----------------------------------------------------------------------
! A=[m], This function is only used to generate the test data.
!   [m+1]
!   [m+2]
! -----------------------------------------------------------------------
function dm_to(m,n) result(A)
    implicit none
    integer,   intent(in)  	::  m,n
    type(Matrix)           	::  A
    integer					::  ierr
	A= dm_seqs(n-m+1,1)+m
    call dm_set_implicit(A,ierr)
end function


! -----------------------------------------------------------------------
! The eyes function is used to generate the simple and complex identity matrixs. 
! For example, if A is a 2*6 matrix, we can use mat_eye(A,ierr) to obtain 
! A= [1 0 1 0 1 0]
!	 [0 1 0 1 1 0]
! if A is a 6*2 matrix, then mat_eye(A,ierr) will generate
! A= [1 0]
!	 [0 1]
!	 [1 0]
!    [0 1]
!	 [1 0]
!    [0 1]
! -----------------------------------------------------------------------
function dm_eyes(m,n) result(A)
	implicit none
	integer,    intent(in)	::	m,n	
	type(Matrix)			::	A
	integer					::	ierr
    
    call mat_eyes(A%x,m,n,ierr)
    call dm_set_implicit(A,ierr)
end function 


! -----------------------------------------------------------------------
! B=A. This function uses the implicit matrix A directly because A is not need to free. 
! -----------------------------------------------------------------------
subroutine dm_copy(B,A)
    implicit none
    type(Matrix),  intent(in)    ::  A
    type(Matrix),  intent(inout) ::  B
    type(Matrix)                 ::  W
    integer						 ::  ierr
	!print *,"B Type=",B%xtype	
	!print *,"A Type=",A%xtype	
	if(B%xtype==MAT_XTYPE_EXPLICIT) then
        W%x=B%x
    endif
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        B%x=A%x
    else
        call mat_copy(A%x,B%x,ierr)
    endif
    if(B%xtype==MAT_XTYPE_EXPLICIT) then
    	!Free the space of B matrix 
        call mat_destroy(W%x,ierr)
    endif
    call dm_set_explicit(B,ierr)
end subroutine


! -----------------------------------------------------------------------
! C=A+B
! -----------------------------------------------------------------------
function dm_add1(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)                ::	C
	integer						::	ierr
	
    call mat_add(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 

function dm_add2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer						::	ierr
	
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    call mat_add(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
   	
	call mat_destroy(B%x,ierr) 
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 

function dm_add3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,			intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_add2(A,real(alpha,8))	
    C%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_add4(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_add2(A,alpha)	
    C%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_add5(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,			intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_add2(A,real(alpha,8))	
end function 

function dm_add6(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,		intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_add2(A,real(alpha,8))	
    C%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_add7(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,		intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_add2(A,real(alpha,8))	
    C%xtype=MAT_XTYPE_IMPLICIT 
end function 

! -----------------------------------------------------------------------
! C=A-B
! -----------------------------------------------------------------------
function dm_minus1(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_minus(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 

function dm_minus2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    call mat_minus(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
   	
	call mat_destroy(B%x,ierr) 
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 

function dm_minus3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,			intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_minus2(A,real(alpha,8))	
	C%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_minus4(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer						::	ierr
	
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    call mat_minus(B%x,A%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
   	
	call mat_destroy(B%x,ierr) 
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 

function dm_minus5(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,			intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_minus4(real(alpha,8),A)	
    C%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_minus6(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,		intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_minus2(A,real(alpha,8))	
    C%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_minus7(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,		intent(in)	::  alpha 
	type(Matrix)                ::	C
	C=dm_minus4(real(alpha,8),A)	
    C%xtype=MAT_XTYPE_IMPLICIT 
end function 

! -----------------------------------------------------------------------
! C=[A B] 
! -----------------------------------------------------------------------
function dm_hjoin(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_hjoin(A%x,B%x,C%x,ierr)
	call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 

! -----------------------------------------------------------------------
! C=[A] 
!   [B] 
! -----------------------------------------------------------------------
function dm_vjoin(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_vjoin(A%x,B%x,C%x,ierr)
	call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 



! -----------------------------------------------------------------------
! C=A*B
! -----------------------------------------------------------------------
function dm_mult1(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_mult(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 

function dm_mult2(alpha,A) result(B)
	implicit none
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_copy(A%x,B%x,ierr) 
    call mat_scale(B%x,alpha,ierr)
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 

function dm_mult3(A,alpha) result(B)
	implicit none
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              	::	B
   	B=dm_mult2(alpha,A)
	B%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_mult4(alpha,A) result(B)
	implicit none
	real,           intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
   	B=dm_mult2(real(alpha,8),A)
	B%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_mult5(A,alpha) result(B)
	implicit none
	real,           intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
   	B=dm_mult2(real(alpha,8),A)
	B%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_mult6(A,alpha) result(B)
	implicit none
	integer,        intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
   	B=dm_mult2(real(alpha,8),A)
	B%xtype=MAT_XTYPE_IMPLICIT 
end function 

function dm_mult7(alpha,A) result(B)
	implicit none
	integer,        intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
   	B=dm_mult2(real(alpha,8),A)
	B%xtype=MAT_XTYPE_IMPLICIT 
end function 

! -----------------------------------------------------------------------
! C=A.*B
! -----------------------------------------------------------------------
function dm_emult(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_emult(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C=A./B
! -----------------------------------------------------------------------
function dm_ediv(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_ediv(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=repmat(A,m,n)
! -----------------------------------------------------------------------
function dm_rep(A,m,n) result(B)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  m,n 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_rep(A%x,m,n,B%x,ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Sum of elements along with the row or column.
! Suppose A=[1,2,3]
!           [4,5,6],
! then mat_sum(A,1,B) will make B=[5,7,9],
!      mat_sum(A,2,B) will make B=[6 ]
!                                 [15]
! -----------------------------------------------------------------------
function dm_sum(A,ndim) result(B)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  ndim 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_sum(A%x,ndim,B%x,ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Compute Y = a*X + Y.
! -----------------------------------------------------------------------

subroutine dm_axpy1(Y,a,X,ierr)
	implicit none
	type(Matrix),	intent(in)	::  X 
	real(kind=8),   intent(in)	::	a
	type(Matrix), intent(inout) ::  Y 
	integer,		intent(out)	::	ierr
    call mat_axpy(Y%x,a,X%x,ierr)
    call dm_set_explicit(Y,ierr)
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end subroutine

subroutine dm_axpy2(Y,a,X,ierr)
	implicit none
	type(Matrix),	intent(in)	::  X 
	real,           intent(in)	::	a
	type(Matrix), intent(inout) ::  Y 
	integer,		intent(out)	::	ierr
    call mat_axpy(Y%x,real(a,kind=8),X%x,ierr)
    call dm_set_explicit(Y,ierr)
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end subroutine

subroutine dm_axpy3(Y,a,X,ierr)
	implicit none
	type(Matrix),	intent(in)		::  X 
	integer,        intent(in)		::	a
	type(Matrix), 	intent(inout) 	::  Y 
	integer,		intent(out)	::	ierr
    call mat_axpy(Y%x,real(a,kind=8),X%x,ierr)
    call dm_set_explicit(Y,ierr)
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end subroutine


! -----------------------------------------------------------------------
! Compute Y = a*Y + X.
! -----------------------------------------------------------------------
subroutine dm_aypx1(Y,a,X,ierr)
	implicit none
	type(Matrix),	intent(in)	::  X 
	real(kind=8),    intent(in)	::	a
	type(Matrix), intent(inout) ::  Y 
	integer,		intent(out)	::	ierr
    call mat_aypx(Y%x,a,X%x,ierr)
    call dm_set_explicit(Y,ierr)
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end subroutine

subroutine dm_aypx2(Y,a,X,ierr)
	implicit none
	type(Matrix),	intent(in)	::  X 
	real        ,   intent(in)	::	a
	type(Matrix), intent(inout) ::  Y 
	integer,		intent(out)	::	ierr
    call mat_aypx(Y%x,real(a,kind=8),X%x,ierr)
    call dm_set_explicit(Y,ierr)
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end subroutine

subroutine dm_aypx3(Y,a,X,ierr)
	implicit none
	type(Matrix),	intent(in)		::  X 
	integer,   		intent(in)		::	a
	type(Matrix), 	intent(inout) 	::  Y 
	integer,		intent(out)	::	ierr
    call mat_aypx(Y%x,real(a,kind=8),X%x,ierr)
    call dm_set_explicit(Y,ierr)
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end subroutine

! -----------------------------------------------------------------------
! B = A^T.
! -----------------------------------------------------------------------
function dm_trans(A) result(B)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_trans(A%x,B%x,ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C = A*B^T
! -----------------------------------------------------------------------
function dm_xyt(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_xyt(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C = A^T*B
! -----------------------------------------------------------------------
function dm_xty(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_xty(A%x,B%x,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=exp(A) 
! -----------------------------------------------------------------------
function dm_exp(A) result(B)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_math(A%x,MAT_MATH_EXP,B%x,ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=log(A) 
! -----------------------------------------------------------------------
function dm_log(A) result(B)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_math(A%x,MAT_MATH_LOG,B%x,ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=(A).^2 
! -----------------------------------------------------------------------
function dm_squ(A) result(B)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_math(A%x,MAT_MATH_SQU,B%x,ierr)
    call dm_set_implicit(B,ierr)
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=(A).^3 
! -----------------------------------------------------------------------
function dm_cube(A) result(B)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_math(A%x,MAT_MATH_CUBE,B%x,ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 



! -----------------------------------------------------------------------
! B=sqrt(A) 
! -----------------------------------------------------------------------
function dm_sqrt(A) result(B)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              	::	B
	integer						::	ierr
    call mat_math(A%x,MAT_MATH_SQRT,B%x,ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Solve Ax=b 
! -----------------------------------------------------------------------
function dm_solve(A,B) result(X)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)  ::	B
	type(Matrix)            	::	X
	integer						::	ierr
    
    call mat_solve(A%x,B%x,X%x,ierr)
    call dm_set_implicit(X,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Load a standard row-cloumn file into a matrix 
! -----------------------------------------------------------------------
subroutine dm_load(filename,A,ierr)
	implicit none
    character(len=*),   intent(in)  ::  filename 
	type(Matrix),		intent(out)	::  A 
	integer,			intent(out)	::	ierr
    
    call mat_load(filename,A%x,ierr)
	call dm_set_explicit(A,ierr)
end subroutine 


! -----------------------------------------------------------------------
! A(m,n)=value 
! -----------------------------------------------------------------------
subroutine dm_setvalue1(A,m,n,value,ierr)
	implicit none
	type(Matrix)	            ::  A 
	integer,	    intent(in)	::	m,n
	real(kind=8),   intent(in)	::	value
	integer,		intent(out)	::	ierr
	
    call mat_setvalue(A%x,m,n,value,ierr)
end subroutine 

subroutine dm_setvalue2(A,m,n,value,ierr)
	implicit none
	type(Matrix)	            ::  A 
	integer,	    intent(in)	::	m,n
	integer,    	intent(in)	::	value
	integer,		intent(out)	::	ierr
	
    call mat_setvalue(A%x,m,n,real(value,8),ierr)
end subroutine 

subroutine dm_setvalue3(A,m,n,value,ierr)
	implicit none
	type(Matrix)	            ::  A 
	integer,	    intent(in)	::	m,n
	real,    		intent(in)	::	value
	integer,		intent(out)	::	ierr
	
    call mat_setvalue(A%x,m,n,real(value,8),ierr)
end subroutine 


! -----------------------------------------------------------------------
! B=A(rows,cols). Get sub matrix.
! -----------------------------------------------------------------------
function dm_submatrix(A,Rows,Cols) result(B)
	implicit none
	type(Matrix),	intent(in)	::	A
	type(Matrix),	intent(in)	::	Rows
	type(Matrix),	intent(in)	::	Cols
	type(Matrix)				::	B
	integer						::	ierr
    
	call mat_submatrix(A%x,Rows%x,Cols%x,B%x,ierr)
   	call dm_set_implicit(B,ierr)
 
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (Rows%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(Rows%x,ierr)
    endif
    if (Cols%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(Cols%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=A(:,i). Get a column from A.
! -----------------------------------------------------------------------
function dm_getcol(A,n) result(B)
	implicit none
	type(Matrix),	intent(in)	::	A
    integer,        intent(in)  ::  n
	type(Matrix)				::	B
	integer						::	ierr
	type(Matrix)				::  Rows,Cols
	
	Rows=dm_seqs(A%nrow,1)	
	Cols=dm_zeros(1,1)+n	
	
	call mat_submatrix(A%x, Rows%x, Cols%x, B%x, ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
	call mat_destroy(Rows%x,ierr)
	call mat_destroy(Cols%x,ierr)
end function 
 

! -----------------------------------------------------------------------
! B=A(m,:). Get a row from A.
! -----------------------------------------------------------------------
function dm_getrow(A,n) result(B)
	implicit none
	type(Matrix),	intent(in)	::	A
    integer,        intent(in)  ::  n
	type(Matrix)				::	B
	integer						::	ierr
   	type(Matrix)				::  Rows,Cols
	
	Rows=dm_zeros(1,1)+n	
	Cols=dm_seqs(A%ncol,1)	
	
	call mat_submatrix(A%x, Rows%x, Cols%x, B%x, ierr)
    call dm_set_implicit(B,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
	call mat_destroy(Rows%x,ierr)
	call mat_destroy(Cols%x,ierr)
end function 


! -----------------------------------------------------------------------
! Set local values in A.
! -----------------------------------------------------------------------
subroutine dm_setvalues1(A,m,idxm,n,idxn,v,ierr)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	:: 	m,n
	integer,		intent(in)	::	idxm(:),idxn(:)
	real(kind=8),	intent(in)	::	v(:)	
	integer,		intent(out)	::	ierr
   	
	call mat_setvalues(A%x,m,idxm,n,idxn,v,ierr) 
    
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end subroutine 

subroutine dm_setvalues2(A,m,idxm,n,idxn,v,ierr)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	:: 	m,n
	integer,		intent(in)	::	idxm(:),idxn(:)
	real,			intent(in)	::	v(:)	
	integer,		intent(out)	::	ierr

	call mat_setvalues(A%x,m,idxm,n,idxn,real(v,8),ierr) 
    
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end subroutine 

subroutine dm_setvalues3(A,m,idxm,n,idxn,v,ierr)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	:: 	m,n
	integer,		intent(in)	::	idxm(:),idxn(:)
	integer,		intent(in)	::	v(:)	
	integer,		intent(out)	::	ierr

	call mat_setvalues(A%x,m,idxm,n,idxn,real(v,8),ierr) 
    
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end subroutine 


! -----------------------------------------------------------------------
! Get local values in A.
! -----------------------------------------------------------------------
subroutine dm_getvalues(A,m,idxm,n,idxn,v,ierr)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	:: 	m,n
	integer,		intent(in)	::	idxm(:),idxn(:)
	real(kind=8),	intent(inout)	::	v(:)	
	integer,		intent(out)	::	ierr
   	
	call mat_getvalues(A%x,m,idxm,n,idxn,v,ierr) 
    
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end subroutine 

! -----------------------------------------------------------------------
! Norm(A)
! -----------------------------------------------------------------------
function dm_norm_1(A) result(res)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
	type(Matrix),	intent(in)	::	A
	real(kind=8) 				::	res	
	integer						::	ierr
	call mat_norm(A%x,NORM_1,res,ierr) 
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 

function dm_norm_2(A) result(res)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
	type(Matrix),	intent(in)	::	A
	real(kind=8) 				::	res	
	integer						::	ierr
	call mat_norm(A%x,NORM_FROBENIUS,res,ierr) 
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 

function dm_norm_inf(A) result(res)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
	type(Matrix),	intent(in)	::	A
	real(kind=8) 				::	res	
	integer						::	ierr
	call mat_norm(A%x,NORM_INFINITY,res,ierr) 
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C= A<B
! -----------------------------------------------------------------------
function dm_lt1(A,B) result(C)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_compare(A%x,B%x,MAT_COMPARE_LT,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
 
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif

end function

function dm_lt2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    C=dm_lt1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_lt3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_lt1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_lt4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_lt1(A,B)
    call dm_set_implicit(C,ierr)
end function 


! -----------------------------------------------------------------------
! C= (A<=B)
! -----------------------------------------------------------------------
function dm_le1(A,B) result(C)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_compare(A%x,B%x,MAT_COMPARE_LE,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
 
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif

end function

function dm_le2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    C=dm_le1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_le3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_le1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_le4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_le1(A,B)
    call dm_set_implicit(C,ierr)
end function 


! -----------------------------------------------------------------------
! C= (A>B)
! -----------------------------------------------------------------------
function dm_gt1(A,B) result(C)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_compare(A%x,B%x,MAT_COMPARE_GT,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
 
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif

end function

function dm_gt2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    C=dm_gt1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_gt3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_gt1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_gt4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_gt1(A,B)
    call dm_set_implicit(C,ierr)
end function 


! -----------------------------------------------------------------------
! C= (A>=B)
! -----------------------------------------------------------------------
function dm_ge1(A,B) result(C)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_compare(A%x,B%x,MAT_COMPARE_GE,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
 
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif

end function

function dm_ge2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    C=dm_ge1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_ge3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_ge1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_ge4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_ge1(A,B)
    call dm_set_implicit(C,ierr)
end function 


! -----------------------------------------------------------------------
! C= (A==B)
! -----------------------------------------------------------------------
function dm_eq1(A,B) result(C)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_compare(A%x,B%x,MAT_COMPARE_EQ,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
 
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif

end function

function dm_eq2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    C=dm_eq1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_eq3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_eq1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_eq4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_eq1(A,B)
    call dm_set_implicit(C,ierr)
end function 


! -----------------------------------------------------------------------
! C= (A/=B)
! -----------------------------------------------------------------------
function dm_nq1(A,B) result(C)
	implicit none
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    call mat_compare(A%x,B%x,MAT_COMPARE_NQ,C%x,ierr)
    call dm_set_implicit(C,ierr)
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
 
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif

end function

function dm_nq2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,alpha,ierr) 	
    C=dm_nq1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_nq3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_nq1(A,B)
    call dm_set_implicit(C,ierr)
end function 

function dm_nq4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
	call mat_ones(B%x,A%nrow,A%ncol,ierr)
	call mat_scale(B%x,real(alpha,kind=8),ierr) 	
    C=dm_nq1(A,B)
    call dm_set_implicit(C,ierr)
end function 


! -----------------------------------------------------------------------
! Cart2sph(A,B)
! -----------------------------------------------------------------------
subroutine dm_cart2sph(A,B,ierr)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(out)	::  B 
	integer,		intent(out)	::	ierr
	call mat_cart2sph(A%x,B%x,ierr)
	call dm_set_explicit(B,ierr)
end subroutine


! -----------------------------------------------------------------------
! Set implicit type and get other information 
! -----------------------------------------------------------------------

subroutine dm_set_implicit(A,ierr)
	implicit none
	type(Matrix),	intent(inout)	::  A 
	integer,		intent(out)		::	ierr
    A%xtype=MAT_XTYPE_IMPLICIT 
	call mat_getsize(A%x,A%nrow,A%ncol,ierr)
	call mat_getownershiprange(A%x,A%ista,A%iend,ierr) 
end subroutine


! -----------------------------------------------------------------------
! Set explicit type and get other information 
! -----------------------------------------------------------------------
subroutine dm_set_explicit(A,ierr)
	implicit none
	type(Matrix),	intent(inout)	::  A 
	integer,		intent(out)		::	ierr
    A%xtype=MAT_XTYPE_EXPLICIT 
	call mat_getsize(A%x,A%nrow,A%ncol,ierr)
	call mat_getownershiprange(A%x,A%ista,A%iend,ierr) 
end subroutine


end module 
