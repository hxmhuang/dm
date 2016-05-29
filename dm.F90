! -----------------------------------------------------------------------
! Distributed Matrix Computing Library
! -----------------------------------------------------------------------
module dm 
    use dm_type
    implicit none

    interface assignment(=)
        module procedure dm_copyEx
        module procedure dm_copyIm 
    end interface

    interface operator (+)
        module procedure dm_add1
        module procedure dm_add2
        module procedure dm_add3
        module procedure dm_add4
    end interface


    interface dm_destroy
        module procedure dm_destroyEx 
        module procedure dm_destroyIm 
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
    PetscBool       ::  input 
    PetscErrorCode  ::  ierr 
    call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,input,PETSC_NULL_BOOL,ierr)
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



! -----------------------------------------------------------------------
!Create a matrix with m*n size
! -----------------------------------------------------------------------
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
    call MatCreate(PETSC_COMM_WORLD,A%x,ierr);
    call MatSetSizes(A%x,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr)
    call MatSetFromOptions(A%x,ierr)
    call MatSetUp(A%x,ierr)
    call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)
end function 

! -----------------------------------------------------------------------
!Destroy a matrix to free the memory
! -----------------------------------------------------------------------
function dm_destroyEx(A) result(ierr)
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
function dm_destroyIm(A) result(ierr)
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
    type(MatrixIm)                ::  A
    PetscErrorCode              ::  ierr
	PetscLogEvent	            ::  ievent

	call PetscLogEventRegister("dm_zeros",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    A=dm_create(m,n)
    call MatZeroEntries(A%x,ierr)
    
    call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr)
end function


function dm_ones(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n
    PetscErrorCode              ::  ierr
    type(MatrixIm)                ::  A

    PetscInt                    ::  ista,iend
    PetscInt,allocatable        ::  idxn(:)
    PetscScalar,allocatable     ::  row(:),results(:)
    integer                     ::  i,j
	PetscLogEvent	            ::  ievent

	call PetscLogEventRegister("dm_ones",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    A=dm_create(m,n)
    call MatGetOwnershipRange(A%x,ista,iend,ierr)
    !print *,">ista=",ista,"iend=",iend,">ncol=",ncol
    allocate(idxn(n),row(n),results(n))

    do i=ista,iend-1
        do j=1,n
            idxn(j)=j-1
            row(j)=1.0
        enddo
        call MatSetValues(A%x,1,i,n,idxn,row,INSERT_VALUES,ierr)
    enddo
    call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)

    deallocate(idxn,row,results)
    
    call PetscLogEventEnd(ievent,ierr)
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
    PetscErrorCode              ::  ierr
    type(MatrixIm)                ::  A

    PetscInt                    ::  ista,iend
    PetscInt,allocatable        ::  idxn(:)
    PetscScalar,allocatable     ::  row(:),results(:)
    integer                     ::  i,j
	PetscLogEvent	            ::  ievent

	call PetscLogEventRegister("dm_seq",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    A=dm_create(m,n)
    call MatGetOwnershipRange(A%x,ista,iend,ierr)
    !print *,">ista=",ista,"iend=",iend,">ncol=",ncol
    allocate(idxn(n),row(n),results(n))

    do i=ista,iend-1
        do j=1,n
            idxn(j)=j-1
			row(j)=i*n+j
        enddo
        call MatSetValues(A%x,1,i,n,idxn,row,INSERT_VALUES,ierr)
    enddo
    call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)

    deallocate(idxn,row,results)
    
    call PetscLogEventEnd(ievent,ierr)
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
	PetscInt					::	nmax, nmin 
	PetscInt					::  ista,iend,xpos,ypos
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:)
	integer 					:: 	i,j,counter
	PetscLogEvent	            ::  ievent

	call PetscLogEventRegister("dm_eye",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    A=dm_create(m,n) 
    nmin=min(m,n)
	nmax=max(m,n)
 	if(mod(nmax,nmin) /= 0) then
		print *, "Error in dm_eye: the size of matrix A should be (NM)*M or M*(NM)"
		stop	
	endif
	call MatGetOwnershipRange(A%x,ista,iend,ierr)
	allocate(idxn(nmax/nmin),row(nmax/nmin))

    do i=ista,iend-1
    	xpos=mod(i,nmin)

        counter=0	
		do j=1,n
    		ypos=mod(j-1,nmin)
    		if(ypos==xpos) then
                counter=counter+1
    			row(counter)=1.0
			    idxn(counter)=j-1
    		    !print *,"i=",i,"j=",j,"xpos=",xpos,"ypos=",ypos,"idxn(",counter,")=",idxn(counter),"row(",j,")=",row(counter) 
    		endif
		enddo
	   	call MatSetValues(A%x,1,i,counter,idxn,row,INSERT_VALUES,ierr)
	enddo
    call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)
	deallocate(idxn,row)
    call PetscLogEventEnd(ievent,ierr)
end function 



! -----------------------------------------------------------------------
! B=A. This function uses the implicit matrix A directly because A is not need to free. 
! -----------------------------------------------------------------------
subroutine dm_copyIm(B,A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(MatrixIm),  intent(in)  ::  A
    type(Matrix),    intent(out) ::  B
    PetscErrorCode               ::  ierr
    !Free the space of B matrix 
    !call MatDestroy(B%x,ierr)
    B%x=A%x
end subroutine


subroutine dm_copyEx(B,A)
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
    call MatDuplicate(A%x,MAT_COPY_VALUES,B%x,ierr)
end subroutine



! -----------------------------------------------------------------------
! C=A+B
! -----------------------------------------------------------------------
subroutine mat_add(A,B,C,ierr) 
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,	       intent(in)	::  A 
	Mat,	       intent(in)	::  B 
	Mat,           intent(out)  ::	C
	PetscErrorCode,intent(out)  ::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
    PetscScalar                 ::  alpha	
    call MatGetSize(A,nrow1,ncol1,ierr)
	call MatGetSize(B,nrow2,ncol2,ierr)
	if(nrow1/=nrow2 .or. ncol1/=ncol2)then
		print *, "Error: matrix A and matrix B should have the same size"
		stop	
	endif
    
    alpha=1.0
    call MatDuplicate(A,MAT_COPY_VALUES,C,ierr)
	call MatAXPY(C,alpha,B,DIFFERENT_NONZERO_PATTERN,ierr)
end subroutine 


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





end module 


