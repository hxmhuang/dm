! -----------------------------------------------------------------------
! Distributed Matrix Wrapper 
! a----------------------------------------------------------------------
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
        module procedure dm_minus8
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
        module procedure dm_nq4
    end interface
  
  	! element multiple
    interface operator (.em.)
        module procedure dm_emult
    end interface
     
    ! element divide 
    interface operator (.ed.)
        module procedure dm_ediv
    end interface
   
    ! join along with x direction 
    interface operator (.xj.)
        module procedure dm_xjoin
    end interface
 
    ! join along with y direction 
    interface operator (.yj.)
        module procedure dm_yjoin
    end interface
  
    ! join along z direction 
    interface operator (.zj.)
        module procedure dm_zjoin
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

    interface dm_setdiag
        module procedure dm_setdiag1
        module procedure dm_setdiag2
        module procedure dm_setdiag3
    end interface
    
	interface assignment(=)
        module procedure dm_copy
    end interface

	integer 		::  dm_nx
	integer 		::	dm_ny
	integer 		::	dm_nz
	type(Matrix)	:: 	DM_ZERO
contains

! -----------------------------------------------------------------------
! Initialize the distributed matrix environment 
! -----------------------------------------------------------------------
subroutine dm_init1(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)  ::  ierr 
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
end subroutine 

subroutine dm_init2(m,n,k,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
	integer,intent(in)	 ::  m,n,k
    integer,intent(out)  ::  ierr 
	dm_nx=m
	dm_ny=n
	dm_nz=k
	call mat_create(DM_ZERO%x,m,n,k,.true.,ierr)
	call mat_zeros(DM_ZERO%x,ierr)
	call mat_assemble(DM_ZERO%x,ierr)
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
	call mat_destroy(DM_ZERO%x,ierr)
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
! Print a matrix on screen
! -----------------------------------------------------------------------
subroutine dm_printf(str,ierr) 
    implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*)            ::  str
    integer,	   intent(out)	::  ierr

    call PetscPrintf(PETSC_COMM_WORLD,str,ierr)   
end subroutine 

! -----------------------------------------------------------------------
! Create a matrix 
! -----------------------------------------------------------------------
subroutine dm_create(A,m,n,k,isGlobal,ierr)
    implicit none
    integer, 	    intent(in)          ::  m,n,k 
    logical,        intent(in),optional ::  isGlobal 
    type(Matrix),   intent(out)         ::  A
    integer,        intent(out)         ::  ierr
     
    if (present(isGlobal)) then
		A%isGlobal=isGlobal
    else
        A%isGlobal=.true.
    endif
	
	if(m==dm_nx .and. n==dm_ny .and. k==dm_nz .and. (A%isGlobal .eqv. .true.)) then 
        call mat_copy(DM_ZERO%x,A%x,ierr)
    else
		call mat_create(A%x,m,n,k,isGlobal,ierr)
	endif	
	A%nx=m
	A%ny=n
	A%nz=k
end subroutine 


! -----------------------------------------------------------------------
! A=0 
! -----------------------------------------------------------------------
function dm_zeros(m,n,k,isGlobal) result(A)
    implicit none
    integer, 	intent(in)          ::  m,n,k 
    logical,    intent(in),optional ::  isGlobal 
    type(Matrix)                    ::  A
    integer					        ::  ierr
   
    if (present(isGlobal)) then
        call dm_create(A,m,n,k,isGlobal,ierr)
   	else
        call dm_create(A,m,n,k,.true.,ierr)
	endif
	call mat_zeros(A%x,ierr)
    call dm_set_implicit(A,ierr)
end function


! -----------------------------------------------------------------------
! A=1. 
! -----------------------------------------------------------------------
function dm_ones(m,n,k,isGlobal) result(A)
    implicit none
    integer, 	intent(in)          ::  m,n,k 
    logical,    intent(in),optional ::  isGlobal 
    type(Matrix)                    ::  A
    integer					        ::  ierr
    
    if (present(isGlobal)) then
        call dm_create(A,m,n,k,isGlobal,ierr)
    else
        call dm_create(A,m,n,k,.true.,ierr)
    endif
    call mat_constants(A%x,A%nx,A%ny,A%nz,real(1.0,8),ierr)
    call dm_set_implicit(A,ierr)
end function


! -----------------------------------------------------------------------
! The eye function is used to generate the simple and complex identity matrixs. 
! For example, if A is a 2*6 matrix, we can use mat_eye(A,ierr) to obtain 
! A= [1 0 0 0 0 0]
!	 [0 1 0 0 0 0]
! if A is a 6*2 matrix, then mat_eye(A,ierr) will generate
! A= [1 0]
!	 [0 1]
!	 [0 0]
!    [0 0]
!	 [0 0]
!    [0 0]
! -----------------------------------------------------------------------
function dm_eye(m,n,k,isGlobal) result(A)
	implicit none
    integer,    intent(in)  	        ::  m,n,k
    logical,    intent(in),optional     ::  isGlobal 
	type(Matrix)			            ::	A
	integer                             ::	nmax, nmin 
	integer					            ::	ierr
    if (present(isGlobal)) then
        call dm_create(A,m,n,k,isGlobal,ierr)
    else
        call dm_create(A,m,n,k,.true.,ierr)
    endif
    call mat_eye(A%x,m,n,k,ierr)
    call dm_set_implicit(A,ierr)
end function 


! -----------------------------------------------------------------------
! B=A. This function uses the implicit matrix A directly because A is not need to release. 
! -----------------------------------------------------------------------
subroutine dm_copy(B,A)
    implicit none
    type(Matrix),  intent(in)    ::  A
    type(Matrix),  intent(inout) ::  B
    type(Matrix)                 ::  W
    integer						 ::  ierr
	
	B%isGlobal=A%isGlobal
	B%nx=A%nx
	B%ny=A%ny
	B%nz=A%nz
	if(B%xtype==MAT_XTYPE_EXPLICIT) then
        W%x=B%x
    endif
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        B%x=A%x
    else
        call mat_copy(A%x,B%x,ierr)
    endif
    if(B%xtype==MAT_XTYPE_EXPLICIT) then
    	!release the original B matrix 
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
	real(kind=8)				::  alpha	
	integer						::	ierr
	alpha=1.0
    
    if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_add: Matrix A and B should have the same distribution.",ierr)
        stop
    endif
    
	if(A%xtype==MAT_XTYPE_IMPLICIT .and. B%xtype==MAT_XTYPE_IMPLICIT) then
		call dm_copy(C,A)
        call mat_axpy(C%x,alpha,B%x,ierr)
        call mat_destroy(B%x,ierr)
	else if(A%xtype==MAT_XTYPE_IMPLICIT .and. B%xtype==MAT_XTYPE_EXPLICIT) then	
		call dm_copy(C,A)
		call mat_axpy(C%x,alpha,B%x,ierr)
	else if(A%xtype==MAT_XTYPE_EXPLICIT .and. B%xtype==MAT_XTYPE_IMPLICIT) then	
		call dm_copy(C,B)
		call mat_axpy(C%x,alpha,A%x,ierr)
	else if(A%xtype==MAT_XTYPE_EXPLICIT .and. B%xtype==MAT_XTYPE_EXPLICIT) then	
		call dm_copy(C,A)
		call mat_axpy(C%x,alpha,B%x,ierr)
	else
		call dm_printf("Error in dm_add: wrong xtype value in this funtion.",ierr)
        stop
	endif 
	
    call dm_set_implicit(C,ierr)
end function 

function dm_add2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
 	call dm_create(C,A%nx,A%ny,A%nz,A%isGlobal,ierr)
    call mat_constants(C%x,C%nx,C%ny,C%nz,alpha,ierr) 
    C=dm_add1(A,C)
    call dm_set_implicit(C,ierr)
end function 

function dm_add3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,			intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_add2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
end function 

function dm_add4(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_add2(A,alpha)	
    call dm_set_implicit(C,ierr)
end function 

function dm_add5(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,			intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_add2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
end function 

function dm_add6(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,		intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_add2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
end function 

function dm_add7(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,		intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_add2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
end function 

! -----------------------------------------------------------------------
! C=A-B
! -----------------------------------------------------------------------
function dm_minus1(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)                ::	C
	real(kind=8)				::  alpha	
	integer						::	ierr
	alpha=-1.0
    
    if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_minus: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

    if(A%xtype==MAT_XTYPE_IMPLICIT .and. B%xtype==MAT_XTYPE_IMPLICIT) then
		call dm_copy(C,A)
		call mat_axpy(C%x,alpha,B%x,ierr)
        call mat_destroy(B%x,ierr)
	else if(A%xtype==MAT_XTYPE_IMPLICIT .and. B%xtype==MAT_XTYPE_EXPLICIT) then	
		call dm_copy(C,A)
		call mat_axpy(C%x,alpha,B%x,ierr)
	else if(A%xtype==MAT_XTYPE_EXPLICIT .and. B%xtype==MAT_XTYPE_IMPLICIT) then	
		call dm_copy(C,B)
		call mat_aypx(C%x,alpha,A%x,ierr)
	else if(A%xtype==MAT_XTYPE_EXPLICIT .and. B%xtype==MAT_XTYPE_EXPLICIT) then	
		call dm_copy(C,A)
		call mat_axpy(C%x,alpha,B%x,ierr)
	else
		print *,"Error in dm_minus: wrong xtype value in this funtion."
	endif 
	call dm_set_implicit(C,ierr)
end function 

function dm_minus2(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr

 	call dm_create(C,A%nx,A%ny,A%nz,A%isGlobal,ierr)
    call mat_constants(C%x,C%nx,C%ny,C%nz,alpha,ierr) 
    C=dm_minus1(A,C)
 	call dm_set_implicit(C,ierr)
end function 

function dm_minus3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,			intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_minus2(A,real(alpha,8))	
	call dm_set_implicit(C,ierr)
end function 

function dm_minus4(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real(kind=8),	intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
   
 	call dm_create(C,A%nx,A%ny,A%nz,A%isGlobal,ierr)
    call mat_constants(C%x,C%nx,C%ny,C%nz,alpha,ierr) 
    C=dm_minus1(C,A)
    call dm_set_implicit(C,ierr)
end function 

function dm_minus5(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,			intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_minus4(real(alpha,8),A)	
	call dm_set_implicit(C,ierr)
end function 

function dm_minus6(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,		intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_minus2(A,real(alpha,8))	
	call dm_set_implicit(C,ierr)
end function 

function dm_minus7(alpha,A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,		intent(in)	::  alpha 
	type(Matrix)                ::	C
	integer						::	ierr
	C=dm_minus4(real(alpha,8),A)	
	call dm_set_implicit(C,ierr)
end function 

function dm_minus8(A) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	C
	integer						::	ierr
	C=(-1.0)*A
	call dm_set_implicit(C,ierr)
end function 


! -----------------------------------------------------------------------
! C=[A B] 
! -----------------------------------------------------------------------
function dm_xjoin(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr

    if((A%nx/=B%nx) .or. (A%nz/=B%nz) .or.  (A%isGlobal .neqv. B%isGlobal)) then
        call dm_printf("Error in dm_xjoin: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

    call mat_xjoin(A%x,A%nx,A%ny,A%nz,B%x,B%nx,B%ny,B%nz,C%x,ierr)
    call dm_set_implicit(C,ierr)
	C%ny=A%ny+B%ny
    C%nx=A%nx
	C%nz=A%nz

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
function dm_yjoin(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
	!print *, "A%ny=",A%ny,"B%ny=",B%ny
	!print *, "A%nx=",A%nx,"B%nx=",B%nx
    if((A%ny/=B%ny) .or. (A%nz/=B%nz) .or.  (A%isGlobal .neqv. B%isGlobal)) then
		print *, "Error in dm_yjoin: Matrix A and Matrix B should have the same distribution."
		stop	
	endif
    
    call mat_yjoin(A%x,A%nx,A%ny,A%nz,B%x,B%nx,B%ny,B%nz,C%x,ierr)
	call dm_set_implicit(C,ierr)
	C%nx=A%nx+B%nx
    C%ny=A%ny
	C%nz=A%nz
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 

! -----------------------------------------------------------------------
! C=[A 0] 
!   [0 B] 
! -----------------------------------------------------------------------
function dm_zjoin(A,B) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              	::	C
	integer						::	ierr
    
	if((A%nx/=B%nx) .or. (A%ny/=B%ny) .or.  (A%isGlobal .neqv. B%isGlobal)) then
		print *, "Error in dm_yjoin: Matrix A and Matrix B should have the same distribution."
		stop	
	endif
    
    call mat_zjoin(A%x,A%nx,A%ny,A%nz,B%x,B%nx,B%ny,B%nz,C%x,ierr)
	call dm_set_implicit(C,ierr)
	C%nx=A%nx
    C%ny=A%ny
	C%nz=A%nz+B%nz
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

    if((A%ny/=B%nx) .or. (A%nz/=B%nz) .or.  (A%isGlobal .neqv. B%isGlobal)) then
 		print *, "Error in dm_mult: the column of A matrix should equal to the row of B matrix, and the number of z dimension should be same."
		stop	
	endif
    
    call mat_mult(A%x,B%x,C%x,ierr)
    C%isGlobal=A%isGlobal
    call dm_set_implicit(C,ierr)
	C%nx=A%nx
	C%ny=B%ny	
	C%nz=A%nz
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
    B%isGlobal=A%isGlobal
	B%nx=A%nx
	B%ny=A%ny
	B%nz=B%nz	
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
	integer						::	ierr
   	B=dm_mult2(alpha,A)
    call dm_set_implicit(B,ierr)
end function 

function dm_mult4(alpha,A) result(B)
	implicit none
	real,           intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
	integer						::	ierr
   	B=dm_mult2(real(alpha,8),A)
    call dm_set_implicit(B,ierr)
end function 

function dm_mult5(A,alpha) result(B)
	implicit none
	real,           intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
	integer						::	ierr
   	B=dm_mult2(real(alpha,8),A)
    call dm_set_implicit(B,ierr)
end function 

function dm_mult6(A,alpha) result(B)
	implicit none
	integer,        intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
	integer						::	ierr
   	B=dm_mult2(real(alpha,8),A)
    call dm_set_implicit(B,ierr)
end function 

function dm_mult7(alpha,A) result(B)
	implicit none
	integer,        intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
	integer						::	ierr
   	B=dm_mult2(real(alpha,8),A)
    call dm_set_implicit(B,ierr)
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
 
    if((A%nx/=B%nx) .or. (A%ny/=B%ny) .or. (A%nz/=B%nz) .or.  (A%isGlobal .neqv. B%isGlobal)) then
 		print *, "Error in dm_emult: the A matrix and B matrix should have the same distribution." 
		stop	
	endif
    
!   if(A%isGlobal .neqv. B%isGlobal) then
!       call dm_printf("Error in dm_emult: Matrix A and B should have the same distribution.",ierr)
!       stop
!   endif

!   if(A%nrow /= B%nrow .or. A%ncol /= B%ncol)then
!   	print *, "Error in dm_emult: Matrix A and Matrix B should have the same size."
!   	stop	
!   endif
!  
!   call dm_create(C,A%nrow,A%ncol,A%isGlobal,ierr) 
!   call mat_emult(A%x,B%x,C%x,ierr)
!   C%isGlobal=A%isGlobal
!   call dm_set_implicit(C,ierr)
!   
!   if (A%xtype==MAT_XTYPE_IMPLICIT) then
!       call mat_destroy(A%x,ierr)
!   endif
!   if (B%xtype==MAT_XTYPE_IMPLICIT) then
!       call mat_destroy(B%x,ierr)
!   endif
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
  
!   if(A%isGlobal .neqv. B%isGlobal) then
!       call dm_printf("Error in dm_emult: Matrix A and B should have the same distribution.",ierr)
!       stop
!   endif

!   if(A%nrow /= B%nrow .or. A%ncol /= B%ncol)then
!   	print *, "Error in dm_emult: Matrix A and Matrix B should have the same size."
!   	stop	
!   endif
!  
!   call dm_create(C,A%nrow,A%ncol,A%isGlobal,ierr) 
!   call mat_ediv(A%x,B%x,C%x,ierr)
!   C%isGlobal=A%isGlobal
!   call dm_set_implicit(C,ierr)
!   
!   if (A%xtype==MAT_XTYPE_IMPLICIT) then
!       call mat_destroy(A%x,ierr)
!   endif
!   if (B%xtype==MAT_XTYPE_IMPLICIT) then
!       call mat_destroy(B%x,ierr)
!   endif
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
    B%isGlobal=A%isGlobal
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
    B%isGlobal=A%isGlobal
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
	!call dm_create(Y,X%nrow,X%ncol,X%isGlobal,ierr)
    call mat_axpy(Y%x,a,X%x,ierr)
    Y%isGlobal=X%isGlobal
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
    Y%isGlobal=X%isGlobal
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
    Y%isGlobal=X%isGlobal
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
    Y%isGlobal=X%isGlobal
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
    Y%isGlobal=X%isGlobal
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
    Y%isGlobal=X%isGlobal
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
    B%isGlobal=A%isGlobal
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
     
    if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_xyt: Matrix A and B should have the same distribution.",ierr)
        stop
    endif
	call mat_xyt(A%x,B%x,C%x,ierr)
    C%isGlobal=A%isGlobal
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
     
    if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_xyt: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

	call mat_xty(A%x,B%x,C%x,ierr)
    C%isGlobal=A%isGlobal
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
    B%isGlobal=A%isGlobal
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
    B%isGlobal=A%isGlobal
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
    B%isGlobal=A%isGlobal
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
    B%isGlobal=A%isGlobal
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
    B%isGlobal=A%isGlobal
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
    
    if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_solve: Matrix A and B should have the same distribution.",ierr)
        stop
    endif
    
	 if(A%nrow /= B%nrow) then
        call dm_printf("Error in dm_solve: Matrix A and B should have the same row size.",ierr)
        stop
    endif
     
	 if(B%ncol /= 1) then
        call dm_printf("Error in dm_solve: the column of B should be 1.",ierr)
        stop
    endif

	call mat_solve(A%x,B%x,X%x,ierr)
	X%isGlobal=A%isGlobal
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
subroutine dm_load(filename,isGlobal,A,ierr)
	implicit none
    character(len=*),   intent(in)  ::  filename 
    logical,   			intent(in)  ::  isGlobal 
	type(Matrix),		intent(out)	::  A 
	integer,			intent(out)	::	ierr
    
	A%isGlobal=isGlobal
    call mat_load(filename,A%isGlobal,A%x,ierr)
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
function dm_getsub(A,rows,cols) result(B)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	::	rows(:),cols(:)
	type(Matrix)				::	B
	integer						::	ierr
    
	!print *, "A%isGLobal=",A%isGlobal
	call mat_getsub(A%x,rows,cols,B%x,ierr)
	B%isGLobal=A%isGlobal
   	call dm_set_implicit(B,ierr)
 
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
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
	integer 					:: 	i
	
	call mat_getsub(A%x, (/(i,i=0,A%nrow-1)/), (/n/), B%x, ierr)
	B%isGlobal=A%isGlobal
    call dm_set_implicit(B,ierr)
    
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
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
	integer,allocatable			::	rows(:),cols(:)
	integer 					:: 	i
	
	allocate(rows(1),cols(A%ncol))	
	rows(1)=n
	do i=1,A%ncol
		cols(i)=i-1
	enddo	
	
	call mat_getsub(A%x, rows, cols, B%x, ierr)
	B%isGlobal=A%isGlobal
    call dm_set_implicit(B,ierr)
   	
	deallocate(rows,cols)	

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Set local values in A.
! -----------------------------------------------------------------------
subroutine dm_setvalues1(A,idxm,idxn,v,ierr)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	::	idxm(:),idxn(:)
	real(kind=8),	intent(in)	::	v(:)	
	integer,		intent(out)	::	ierr
   	
	call mat_setvalues(A%x,idxm,idxn,v,ierr) 
    
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end subroutine 

subroutine dm_setvalues2(A,idxm,idxn,v,ierr)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	::	idxm(:),idxn(:)
	real,			intent(in)	::	v(:)	
	integer,		intent(out)	::	ierr

	call mat_setvalues(A%x,idxm,idxn,real(v,8),ierr) 
    
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end subroutine 

subroutine dm_setvalues3(A,idxm,idxn,v,ierr)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	::	idxm(:),idxn(:)
	integer,		intent(in)	::	v(:)	
	integer,		intent(out)	::	ierr

	call mat_setvalues(A%x,idxm,idxn,real(v,8),ierr) 
    
	if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end subroutine 


! -----------------------------------------------------------------------
! Get local values in A.
! -----------------------------------------------------------------------
subroutine dm_getvalues(A,idxm,idxn,v,ierr)
	implicit none
	type(Matrix),	intent(in)	::	A
	integer,		intent(in)	::	idxm(:)
	integer,		intent(in)	::	idxn(:)
	real(kind=8),	intent(inout)	::	v(:)	
	integer,		intent(out)	::	ierr
   	
	call mat_getvalues(A%x,idxm,idxn,v,ierr) 
    
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
	
	if(A%nrow/=B%nrow .or. A%ncol/=B%ncol)then
		print *, "Error: Matrix A and matrix B should have the same size"
		stop	
	endif
 
	if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_lt: Matrix A and B should have the same distribution.",ierr)
        stop
    endif
    
	call mat_compare(A%x,B%x,MAT_COMPARE_LT,C%x,ierr)
	C%isGlobal=A%isGlobal
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
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,alpha,ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_lt1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_lt3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_lt1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_lt4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	 
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!
!!!!C=dm_lt1(A,B)
!!!!call dm_set_implicit(C,ierr)
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
 
	if(A%nrow/=B%nrow .or. A%ncol/=B%ncol)then
		print *, "Error in dm_le: Matrix A and matrix B should have the same size"
		stop	
	endif
 
	if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_le: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_LE,C%x,ierr)
	C%isGlobal=A%isGlobal	
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
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,alpha,ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_le1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_le3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_le1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_le4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_le1(A,B)
!!!!call dm_set_implicit(C,ierr)
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
	
	if(A%nrow/=B%nrow .or. A%ncol/=B%ncol)then
		print *, "Error in dm_gt: Matrix A and matrix B should have the same size"
		stop	
	endif
 
	if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_gt: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_GT,C%x,ierr)
	C%isGlobal=A%isGlobal	
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
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,alpha,ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)

!!!!C=dm_gt1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_gt3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_gt1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_gt4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_gt1(A,B)
!!!!call dm_set_implicit(C,ierr)
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
	
	if(A%nrow/=B%nrow .or. A%ncol/=B%ncol)then
		print *, "Error in dm_ge: Matrix A and matrix B should have the same size"
		stop	
	endif
 
	if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_ge: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_GE,C%x,ierr)
	C%isGlobal=A%isGlobal	
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
	integer						::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,alpha,ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_ge1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_ge3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer						::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_ge1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_ge4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer						::	ierr
!	call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!   call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!   B%isGlobal=A%isGlobal
!   call dm_set_implicit(B,ierr)
!   C=dm_ge1(A,B)
!   call dm_set_implicit(C,ierr)
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
	
	if(A%nrow/=B%nrow .or. A%ncol/=B%ncol)then
		print *, "Error in dm_eq: Matrix A and matrix B should have the same size"
		stop	
	endif
 
	if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_eq: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_EQ,C%x,ierr)
	C%isGlobal=A%isGlobal
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
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,alpha,ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_eq1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_eq3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_eq1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 

function dm_eq4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_eq1(A,B)
!!!!call dm_set_implicit(C,ierr)
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
	
	if(A%nrow/=B%nrow .or. A%ncol/=B%ncol)then
		print *, "Error in dm_nq: Matrix A and matrix B should have the same size"
		stop	
	endif
 
	if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_nq: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_NQ,C%x,ierr)
	C%isGlobal=A%isGlobal
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
!	call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!   call mat_constants(B%x,alpha,ierr) 	
!   B%isGlobal=A%isGlobal
!   call dm_set_implicit(B,ierr)
!   C=dm_nq1(A,B)
!   call dm_set_implicit(C,ierr)
end function 

function dm_nq3(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	real,	        intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!	call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!   call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!   B%isGlobal=A%isGlobal
!   call dm_set_implicit(B,ierr)
!   C=dm_nq1(A,B)
!   call dm_set_implicit(C,ierr)
end function 

function dm_nq4(A,alpha) result(C)
	implicit none
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  alpha 
	type(Matrix)                ::	B
	type(Matrix)                ::	C
	integer::	ierr
!!!!call mat_create(B%x,A%nrow,A%ncol,A%isGlobal,ierr) 	
!!!!call mat_constants(B%x,real(alpha,kind=8),ierr) 	
!!!!B%isGlobal=A%isGlobal
!!!!call dm_set_implicit(B,ierr)
!!!!C=dm_nq1(A,B)
!!!!call dm_set_implicit(C,ierr)
end function 


! -----------------------------------------------------------------------
! sparse(i,j,A,B,ierr)
! -----------------------------------------------------------------------
function dm_sparse(Ind_i,Ind_j,A,m,n) result(B)
	implicit none
	type(Matrix),	intent(in)	::  Ind_i,Ind_j 
	type(Matrix),	intent(in)	::  A 
	integer,	intent(in)		::	m,n	
	type(Matrix)				::  B 
	integer						::	ierr
 
!!!!if((Ind_i%isGlobal .neqv. Ind_j%isGlobal) .or. (A%isGlobal .neqv. Ind_i%isGlobal) ) then
!!!!    call dm_printf("Error in dm_sparse: Matrix Ind_i, Ind_j and A should have the same distribution.",ierr)
!!!!    stop
!!!!endif
!!!!B=dm_zeros(m,n,A%isGlobal)
!!!!call mat_sparse(Ind_i%x,Ind_j%x,A%x,m,n,B%x,ierr)
!!!!call dm_set_implicit(B,ierr)
end function 



! -----------------------------------------------------------------------
! Cart2sph(A,B)
! -----------------------------------------------------------------------
subroutine dm_cart2sph(A,B,ierr)
	implicit none
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(out)	::  B 
	integer,		intent(out)	::	ierr

!   B=dm_zeros(A%nrow,A%ncol,A%isGlobal)
!   call mat_cart2sph(A%x,B%x,ierr)
!   call dm_set_explicit(B,ierr)
end subroutine


! -----------------------------------------------------------------------
! Set implicit type and get other information 
! -----------------------------------------------------------------------

subroutine dm_set_implicit(A,ierr)
	implicit none
	type(Matrix),	intent(inout)	::  A 
	integer,		intent(out)		::	ierr
    A%xtype=MAT_XTYPE_IMPLICIT 
	!call mat_getsize(A%x,A%nrow,A%ncol,ierr)
	!call mat_getownershiprange(A%x,A%ista,A%iend,ierr) 
end subroutine


! -----------------------------------------------------------------------
! Set explicit type and get other information 
! -----------------------------------------------------------------------
subroutine dm_set_explicit(A,ierr)
	implicit none
	type(Matrix),	intent(inout)	::  A 
	integer,		intent(out)		::	ierr
    A%xtype=MAT_XTYPE_EXPLICIT 
	!call mat_getsize(A%x,A%nrow,A%ncol,ierr)
	!call mat_getownershiprange(A%x,A%ista,A%iend,ierr) 
end subroutine


! -----------------------------------------------------------------------
! Set the diagonal of A to constant value 
! -----------------------------------------------------------------------
subroutine dm_setdiag1(A,value,ierr) 
	implicit none
	type(Matrix),	intent(inout)	::  A 
	real(kind=8),	intent(in)		::  value
	integer							::	ierr
	
	call mat_setdiag(A%x,value,ierr)
    call dm_set_explicit(A,ierr)
end subroutine

subroutine dm_setdiag2(A,value,ierr) 
	implicit none
	type(Matrix),	intent(inout)	::  A 
	real,			intent(in)		::  value
	integer							::	ierr
	
	call mat_setdiag(A%x,real(value,kind=8),ierr)
    call dm_set_explicit(A,ierr)
end subroutine

subroutine dm_setdiag3(A,value,ierr) 
	implicit none
	type(Matrix),	intent(inout)	::  A 
	integer,		intent(in)		::  value
	integer							::	ierr
	
	call mat_setdiag(A%x,real(value,kind=8),ierr)
    call dm_set_explicit(A,ierr)
end subroutine


! -----------------------------------------------------------------------
! Set one column of A to another matrix B 
! -----------------------------------------------------------------------
subroutine dm_setcol(A,idxn,B,ierr)
	implicit none

	type(Matrix),	intent(inout)	::  A
	integer,		intent(in)		:: 	idxn 
	type(Matrix),	intent(in)		::  B 
	integer,		intent(out)		::	ierr

    if(A%isGlobal .neqv. B%isGlobal) then
        call dm_printf("Error in dm_setcol: Matrix A and B should have the same distribution.",ierr)
        stop
    endif

	if(A%nrow/= B%nrow)then
		print *, "Error in dm_setcol: Matrix A and Matrix B should have the same row size."
		stop	
	endif

	call mat_setcol(A%x,idxn,B%x,ierr)	
    call dm_set_explicit(A,ierr)
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end subroutine


subroutine dm_test(m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
 	PetscInt,		intent(in)	::	m,n	
    PetscErrorCode,	intent(out)	::	ierr
!   Mat				            ::	A
!   Vec							:: 	b,x
!   KSP                         ::  ksp
!   PC                          ::  pc
!	PetscInt					::  ista1,iend1
!	PetscInt,allocatable		::	idxn(:)
!	PetscScalar,allocatable		::	row(:)
!	integer 					:: 	i,j
!	PetscScalar					:: 	alpha
end subroutine

end module 
