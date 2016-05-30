! -----------------------------------------------------------------------
! Distributed Matrix Computing Interface 
! -----------------------------------------------------------------------
module dm 
    use dm_type
    use dm_mat
    implicit none

    interface assignment(=)
        module procedure dm_copy1
        module procedure dm_copy2
    end interface

    interface operator (+)
        module procedure dm_add1
        module procedure dm_add2
        module procedure dm_add3
        module procedure dm_add4
    end interface

    interface operator (*)
        module procedure dm_mult1
        module procedure dm_mult2
        module procedure dm_mult3
        module procedure dm_mult4
    end interface

    interface operator (.eprod.)
        module procedure dm_eprod1
        module procedure dm_eprod2
        module procedure dm_eprod3
        module procedure dm_eprod4
    end interface

    interface operator (.hjoin.)
        module procedure dm_hjoin1
        module procedure dm_hjoin2
        module procedure dm_hjoin3
        module procedure dm_hjoin4
    end interface

    interface dm_rep
        module procedure dm_rep1
        module procedure dm_rep2
    end interface

    interface dm_sum
        module procedure dm_sum1
        module procedure dm_sum2
    end interface


    interface dm_destroy
        module procedure dm_destroy1 
        module procedure dm_destroy2 
    end interface

contains

! -----------------------------------------------------------------------
! Initialize the distributed matrix environment 
! -----------------------------------------------------------------------
function dm_init() result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscErrorCode  ::  ierr 
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
end function

! -----------------------------------------------------------------------
! Get the rank number of the current process in the commmunicator 
! -----------------------------------------------------------------------
function dm_comm_rank() result(myrank)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer         ::  myrank
    PetscErrorCode  ::  ierr 
	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
end function


! -----------------------------------------------------------------------
! Get the size of processes in the commmunicator 
! -----------------------------------------------------------------------
function dm_comm_size() result(mysize)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer         ::  mysize
    PetscErrorCode  ::  ierr 
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)
end function


! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
function dm_get_int(str) result(input)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*)::  str
    PetscInt        ::  input 
    PetscErrorCode  ::  ierr 
    call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,input,PETSC_NULL_BOOL,ierr)
end function

! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
function dm_get_bool(str) result(input)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*)::  str
    logical         ::  input 
    PetscBool       ::  flag 
    PetscErrorCode  ::  ierr
    call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,flag,PETSC_NULL_BOOL,ierr)
    print *,"str=",str," flag=",flag 
    if(flag .eqv. PETSC_TRUE) then
        input=.true.
    else
        input=.false.
    endif
end function


! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
function dm_get_real(str) result(input)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*)::  str
    PetscReal       ::  input 
    PetscErrorCode  ::  ierr 
    call PetscOptionsGetReal(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,input,PETSC_NULL_BOOL,ierr)
end function


! -----------------------------------------------------------------------
! Finalize the distributed matrix environment 
! -----------------------------------------------------------------------
function dm_finalize() result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscErrorCode  ::  ierr 
    call PetscFinalize(ierr)
end function


function dm_create(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n 
    type(MatrixIm)              ::  A
    PetscErrorCode              ::  ierr
    ! generate matrix A with size m*n
    call mat_create(A%x,m,n,ierr)
end function 

! -----------------------------------------------------------------------
!Destroy a matrix to free the memory
! -----------------------------------------------------------------------
function dm_destroy1(A) result(ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),   intent(in)  ::  A
    PetscErrorCode              ::  ierr
    ! destroy matrix A
    call MatDestroy(A%x,ierr)
end function 

! -----------------------------------------------------------------------
!Destroy a implict matrix to free the memory
! -----------------------------------------------------------------------
function dm_destroy2(A) result(ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(MatrixIm),   intent(in)  ::  A
    PetscErrorCode              ::  ierr
    ! destroy matrix A
    call MatDestroy(A%x,ierr)
end function 

! -----------------------------------------------------------------------
! Print a matrix on screen
! -----------------------------------------------------------------------
function dm_view(A) result(ierr) 
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),  intent(in)   ::  A
    PetscErrorCode              ::  ierr
	call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatView(A%x,PETSC_VIEWER_STDOUT_WORLD, ierr)
end function 

! -----------------------------------------------------------------------
! A=0 
! -----------------------------------------------------------------------
function dm_zeros(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n 
    type(MatrixIm)              ::  A
    PetscErrorCode              ::  ierr
    
    call mat_zeros(A%x,m,n,ierr)
end function



function dm_ones(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n
    type(MatrixIm)              ::  A
    PetscErrorCode              ::  ierr
    
    call mat_ones(A%x,m,n,ierr)
end function


! -----------------------------------------------------------------------
! A=[1 2 3], This function is only used to generate the test data.
!   [4 5 6]
!   [7 8 9]
! -----------------------------------------------------------------------
function dm_seqs(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n
    type(MatrixIm)              ::  A
    PetscErrorCode              ::  ierr

    call mat_seqs(A%x,m,n,ierr)
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
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	PetscInt,       intent(in)	::	m,n	
	type(MatrixIm)			    ::	A
	PetscErrorCode	            ::	ierr
    
    call mat_eyes(A%x,m,n,ierr)
end function 



! -----------------------------------------------------------------------
! B=A. This function uses the implicit matrix A directly because A is not need to free. 
! -----------------------------------------------------------------------
subroutine dm_copy1(B,A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(MatrixIm),  intent(in)  ::  A
    type(Matrix),    intent(out) ::  B
    !PetscErrorCode               ::  ierr
    !Free the space of B matrix 
    !call MatDestroy(B%x,ierr)
    B%x=A%x
end subroutine


subroutine dm_copy2(B,A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),   intent(in)  ::  A
    type(Matrix),   intent(out) ::  B
    PetscErrorCode              ::  ierr
    !Free the space of B matrix 
    !call MatDestroy(B%x,ierr)
    call mat_copy(A%x,B%x,ierr)
end subroutine



! -----------------------------------------------------------------------
! C=A+B
! -----------------------------------------------------------------------
function dm_add1(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_add(A%x,B%x,C%x,ierr)
end function 


function dm_add2(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	type(MatrixIm),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_add(A%x,B%x,C%x,ierr)
    call MatDestroy(A%x,ierr)
    call MatDestroy(B%x,ierr)
end function 


function dm_add3(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_add(A%x,B%x,C%x,ierr)
    call MatDestroy(A%x,ierr)
end function 


function dm_add4(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(MatrixIm),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_add(A%x,B%x,C%x,ierr)
    call MatDestroy(B%x,ierr)
end function 


! -----------------------------------------------------------------------
! C=[A B] 
! -----------------------------------------------------------------------
function dm_hjoin1(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_hjoin(A%x,B%x,C%x,ierr)
end function 


function dm_hjoin2(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	type(MatrixIm),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_hjoin(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(A)
    ierr=dm_destroy(B)
end function 


function dm_hjoin3(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_hjoin(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(A)
end function 


function dm_hjoin4(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(MatrixIm),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_hjoin(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(B)
end function 


! -----------------------------------------------------------------------
! C=A*B
! -----------------------------------------------------------------------
function dm_mult1(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_mult(A%x,B%x,C%x,ierr)
end function 


function dm_mult2(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	type(MatrixIm),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_mult(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(A)
    ierr=dm_destroy(B)
end function 


function dm_mult3(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_mult(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(A)
end function 


function dm_mult4(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(MatrixIm),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_mult(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(B)
end function 


! -----------------------------------------------------------------------
! B=A1.*A2
! -----------------------------------------------------------------------
function dm_eprod1(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_eprod(A%x,B%x,C%x,ierr)
end function 


function dm_eprod2(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	type(MatrixIm),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_eprod(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(A)
    ierr=dm_destroy(B)
end function 


function dm_eprod3(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_eprod(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(A)
end function 


function dm_eprod4(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(MatrixIm),	intent(in)	::  B 
	type(MatrixIm)              ::	C
	PetscErrorCode      		::	ierr
    call mat_eprod(A%x,B%x,C%x,ierr)
    ierr=dm_destroy(B)
end function 


! -----------------------------------------------------------------------
! B=repmat(A,m,n)
! -----------------------------------------------------------------------
function dm_rep1(A,m,n) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	integer,	    intent(in)	::  m,n 
	type(MatrixIm)              ::	B
	PetscErrorCode      		::	ierr
    call mat_rep(A%x,m,n,B%x,ierr)
    ierr=dm_destroy(A)
end function 


function dm_rep2(A,m,n) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  m,n 
	type(MatrixIm)              ::	B
	PetscErrorCode      		::	ierr
    call mat_rep(A%x,m,n,B%x,ierr)
end function 


! -----------------------------------------------------------------------
! Sum of elements along with the row or column.
! Suppose A=[1,2,3]
!           [4,5,6],
! then mat_sum(A,1,B) will make B=[5,7,9],
!      mat_sum(A,2,B) will make B=[6 ]
!                                 [15]
! -----------------------------------------------------------------------
function dm_sum1(A,ndim) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  ndim 
	type(MatrixIm)              ::	B
	PetscErrorCode      		::	ierr
    call mat_sum(A%x,ndim,B%x,ierr)
end function 


function dm_sum2(A,ndim) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(MatrixIm),	intent(in)	::  A 
	integer,	    intent(in)	::  ndim 
	type(MatrixIm)              ::	B
	PetscErrorCode      		::	ierr
    call mat_sum(A%x,ndim,B%x,ierr)
    ierr=dm_destroy(A)
end function 



end module 


