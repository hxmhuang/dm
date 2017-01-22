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
  end interface operator (+)

  interface operator (-)
     module procedure dm_minus1
     module procedure dm_minus2
     module procedure dm_minus3
     module procedure dm_minus4
     module procedure dm_minus5
     module procedure dm_minus6
     module procedure dm_minus7
     module procedure dm_minus8
  end interface operator (-)

  interface operator (*)
     module procedure dm_mult1
     module procedure dm_mult2
     module procedure dm_mult3
     module procedure dm_mult4
     module procedure dm_mult5
     module procedure dm_mult6
     module procedure dm_mult7
  end interface operator (*)

  interface operator (<)
     module procedure dm_lt1
     module procedure dm_lt2
     module procedure dm_lt3
     module procedure dm_lt4
  end interface operator (<)

  interface operator (<=)
     module procedure dm_le1
     module procedure dm_le2
     module procedure dm_le3
     module procedure dm_le4
  end interface operator (<=)

  interface operator (>)
     module procedure dm_gt1
     module procedure dm_gt2
     module procedure dm_gt3
     module procedure dm_gt4
  end interface operator (>)

  interface operator (>=)
     module procedure dm_ge1
     module procedure dm_ge2
     module procedure dm_ge3
     module procedure dm_ge4
  end interface operator (>=)

  interface operator (==)
     module procedure dm_eq1
     module procedure dm_eq2
     module procedure dm_eq3
     module procedure dm_eq4
  end interface operator (==)

  interface operator (/=)
     module procedure dm_nq1
     module procedure dm_nq2
     module procedure dm_nq3
     module procedure dm_nq4
  end interface operator (/=)

  ! element multiple
  interface operator (.em.)
     module procedure dm_emult
  end interface operator (.em.)

  ! element divide 
  interface operator (.ed.)
     module procedure dm_ediv
  end interface operator (.ed.)

  ! join along with x direction 
  interface operator (.xj.)
     module procedure dm_xjoin
  end interface operator (.xj.)

  ! join along with y direction 
  interface operator (.yj.)
     module procedure dm_yjoin
  end interface operator (.yj.)

  ! join along z direction 
  interface operator (.zj.)
     module procedure dm_zjoin
  end interface operator (.zj.)

  ! INV(A)*B
  interface operator (.inv.)
     module procedure dm_solve
  end interface operator (.inv.)

  interface dm_axpy
     module procedure dm_axpy1
     module procedure dm_axpy2
     module procedure dm_axpy3
  end interface dm_axpy

  interface dm_aypx
     module procedure dm_aypx1
     module procedure dm_aypx2
     module procedure dm_aypx3
  end interface dm_aypx

  interface dm_setvalue
     module procedure dm_setvalue1
     module procedure dm_setvalue2
     module procedure dm_setvalue3
  end interface dm_setvalue

  interface dm_setvalues
     module procedure dm_setvalues1
     module procedure dm_setvalues2
     module procedure dm_setvalues3
  end interface dm_setvalues

  interface dm_setdiag
     module procedure dm_setdiag1
     module procedure dm_setdiag2
     module procedure dm_setdiag3
  end interface dm_setdiag

  interface dm_gendiag
     module procedure dm_gendiag1
     module procedure dm_gendiag2
     module procedure dm_gendiag3
  end interface dm_gendiag
  
  interface dm_pow
     module procedure dm_pow1
     module procedure dm_pow2
     module procedure dm_pow3
  end interface dm_pow
  
  interface assignment(=)
     module procedure dm_copy
  end interface assignment(=)

  interface
     subroutine abort() bind(C, name="abort")
     end subroutine abort
  end interface

  
  integer 	::  dm_nx
  integer 	::  dm_ny
  integer 	::  dm_nz
  type(Matrix)	::  DM_ZERO
contains

  ! -----------------------------------------------------------------------
  ! Initialize the distributed matrix environment 
  ! -----------------------------------------------------------------------
  subroutine dm_init1(ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)  ::  ierr 
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  end subroutine dm_init1

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
  end subroutine dm_init2


  ! -----------------------------------------------------------------------
  ! Get the rank number of the current process in the commmunicator 
  ! -----------------------------------------------------------------------
  subroutine dm_comm_rank(myrank,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)         ::  myrank
    integer,intent(out)  		::  ierr 
    call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
  end subroutine dm_comm_rank


  ! -----------------------------------------------------------------------
  ! Get the size of processes in the commmunicator 
  ! -----------------------------------------------------------------------
  subroutine dm_comm_size(mysize,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)         ::  mysize
    integer,intent(out)  		::  ierr 
    call MPI_Comm_size(PETSC_COMM_WORLD,mysize,ierr)
  end subroutine dm_comm_size


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
  end subroutine dm_option_int


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
  end subroutine dm_option_bool


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
  end subroutine dm_option_real


  ! -----------------------------------------------------------------------
  ! Finalize the distributed matrix environment 
  ! -----------------------------------------------------------------------
  subroutine dm_finalize(ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)	::  ierr 
    call mat_destroy(DM_ZERO%x,ierr)
    call PetscFinalize(ierr)
  end subroutine dm_finalize

  subroutine dm_abort(ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
    integer,intent(out)	 ::  ierr 

    call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
    
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

  ! subroutine dm_destroys(A,ierr)
  !   implicit none
  !   type(Matrix),   intent(in)  ::  A(:)
  !   integer,intent(out)			::  ierr
  !   integer :: i
    
  !   do i=1,size(A)
  !      print*, "A(i)",i,A(i)%nx,A(i)%ny,A(i)%nz,A(i)%x
  !      !call mat_destroy(A(i)%x,ierr)       
  !   enddo
  ! end subroutine 

  ! -----------------------------------------------------------------------
  ! Print a matrix on screen
  ! -----------------------------------------------------------------------
  subroutine dm_view(A,ierr) 
    implicit none
    type(Matrix),  intent(in)   ::  A
    integer,	   intent(out)	::  ierr 
    call mat_view(A%x,ierr)

    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif
  end subroutine dm_view


  ! -----------------------------------------------------------------------
  ! Print a matrix on screen
  ! -----------------------------------------------------------------------
  subroutine dm_printf(str,ierr) 
    implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*)            ::  str
    integer,	   intent(out)	::  ierr

    call PetscPrintf(PETSC_COMM_WORLD,str,ierr)   
  end subroutine dm_printf

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
  end subroutine dm_create

  subroutine dm_print_info(A, ierr)
    implicit none
    type(Matrix), intent(in) :: A
    integer, intent(out) :: ierr
    integer :: myrank

    call dm_comm_rank(myrank, ierr)

    if(myrank == 0) then
       write(*, '(A,A,I3,A,I3,A,I3,A,L3)') &
            "Matrix Info:", "  nx=", A%nx, "  ny=", A%ny, "  nz=", A%nz, &
            "  isGlobal=", A%isGlobal
    endif
    
  end subroutine dm_print_info
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
  end function dm_zeros


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
  end function dm_ones

  function dm_constants(m,n,k, value, isGlobal) result(A)
    implicit none
    integer,    intent(in)  	        ::  m,n,k
    logical,    intent(in),optional     ::  isGlobal 
    type(Matrix)			            ::	A
    integer					            ::	ierr
    real(kind=8) :: value
    if (present(isGlobal)) then
       call dm_create(A,m,n,k,isGlobal,ierr)
    else
       call dm_create(A,m,n,k,.true.,ierr)
    endif
    call mat_constants(A%x,m,n,k,value,ierr)
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
    integer					            ::	ierr
    if (present(isGlobal)) then
       call dm_create(A,m,n,k,isGlobal,ierr)
    else
       call dm_create(A,m,n,k,.true.,ierr)
    endif
    call mat_eye(A%x,m,n,k,ierr)
    call dm_set_implicit(A,ierr)
  end function 

  !-------------------------------------------------------------------------
  !> Generate a number sequence start from 0 with a dimension of (m,n,k)
  !-------------------------------------------------------------------------
  function dm_seqs(m,n,k,isGlobal) result(A)
    implicit none
    integer,    intent(in)  	        ::  m,n,k
    logical,    intent(in),optional     ::  isGlobal 
    type(Matrix)			            ::	A
    integer					            ::	ierr
    if (present(isGlobal)) then
       call dm_create(A,m,n,k,isGlobal,ierr)
    else
       call dm_create(A,m,n,k,.true.,ierr)
    endif
    call mat_seqs(A%x,m,n,k,ierr)
    call dm_set_implicit(A,ierr)
  end function 

  !-------------------------------------------------------------------------
  !> Generate a random number sequence ranging (0,1) with a dimension of (m,n,k)
  !-------------------------------------------------------------------------
  function dm_rand(m,n,k,isGlobal) result(A)
    implicit none
    integer,    intent(in)  	        ::  m,n,k
    logical,    intent(in),optional     ::  isGlobal 
    type(Matrix)	                ::	A
    integer					            ::	ierr
    if (present(isGlobal)) then
       call dm_create(A,m,n,k,isGlobal,ierr)
    else
       call dm_create(A,m,n,k,.true.,ierr)
    endif
    call mat_rand(A%x,m,n,k,ierr)
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

    if(B%nx /= 0 .and. B%nx /= A%nx) then
       print*, "WARNING: Reassign matrix with different size."
    endif
    
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
  end subroutine dm_copy


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
       call dm_printf("Error in dm_add: Matrix A and B &
            &should have the same distribution.",ierr)
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
  end function dm_add1

  function dm_add2(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer		::	ierr
    PetscScalar :: val
    call dm_create(C,A%nx,A%ny,A%nz,A%isGlobal,ierr)
    !call mat_constants(C%x,C%nx,C%ny,C%nz,alpha,ierr) 
    !C=dm_add1(A,C)
    
    val = alpha
    
    call mat_apx(C%x, val, A%x, A%nx, A%ny, A%nz, ierr)
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif
    
    call dm_set_implicit(C,ierr)
  end function dm_add2

  function dm_add3(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,			intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_add2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
  end function dm_add3

  function dm_add4(alpha,A) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_add2(A,alpha)	
    call dm_set_implicit(C,ierr)
  end function dm_add4

  function dm_add5(alpha,A) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,			intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_add2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
  end function dm_add5

  function dm_add6(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,		intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_add2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
  end function dm_add6

  function dm_add7(alpha,A) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,		intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_add2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
  end function dm_add7

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
       call dm_printf("Error in dm_minus: Matrix A and B &
            &should have the same distribution.",ierr)
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
  end function dm_minus1

  function dm_minus2(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr

    call dm_create(C,A%nx,A%ny,A%nz,A%isGlobal,ierr)
    !call mat_constants(C%x,C%nx,C%ny,C%nz,alpha,ierr)
    call mat_axpb(C%x, real(1, kind=8), A%x, -alpha, A%nx, A%ny, A%nz, ierr)
    !C=dm_minus1(A,C)
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif
    call dm_set_implicit(C,ierr)
  end function dm_minus2

  function dm_minus3(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,			intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_minus2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
  end function dm_minus3

  function dm_minus4(alpha,A) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr

    call dm_create(C,A%nx,A%ny,A%nz,A%isGlobal,ierr)
    !call mat_constants(C%x,C%nx,C%ny,C%nz,alpha,ierr) 
    !C=dm_minus1(C,A)
    call mat_axpb(C%x, real(-1,kind=8), A%x, alpha, A%nx, A%ny, A%nz, ierr)
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif
    call dm_set_implicit(C,ierr)
  end function dm_minus4

  function dm_minus5(alpha,A) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,			intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_minus4(real(alpha,8),A)	
    call dm_set_implicit(C,ierr)
  end function dm_minus5

  function dm_minus6(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,		intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_minus2(A,real(alpha,8))	
    call dm_set_implicit(C,ierr)
  end function dm_minus6

  function dm_minus7(alpha,A) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,		intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer						::	ierr
    C=dm_minus4(real(alpha,8),A)	
    call dm_set_implicit(C,ierr)
  end function dm_minus7

  function dm_minus8(A) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    type(Matrix)                ::	C
    integer						::	ierr
    C=(-1.0)*A
    call dm_set_implicit(C,ierr)
  end function dm_minus8


  ! -----------------------------------------------------------------------
  ! C=[A] 
  !   [B] 
  ! -----------------------------------------------------------------------
  function dm_xjoin(A,B) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer						::	ierr
    !print *, "A%ny=",A%ny,"B%ny=",B%ny
    !print *, "A%nx=",A%nx,"B%nx=",B%nx
    if((A%ny/=B%ny) .or. (A%ny/=B%ny) .or.  (A%isGlobal .neqv. B%isGlobal)) then
       print *, "Error in dm_xjoin: Matrix A and Matrix B &
            &should have the same distribution."
       stop	
    endif

    call mat_xjoin(A%x,A%nx,A%ny,A%nz,B%x,B%nx,B%ny,B%nz,C%x,ierr)
    call dm_set_implicit(C,ierr)
    C%isGlobal = A%isGlobal
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
  ! C=[A B] 
  ! -----------------------------------------------------------------------
  function dm_yjoin(A,B) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer						::	ierr

    if((A%nx/=B%nx) .or. (A%nz/=B%nz) .or.  (A%isGlobal .neqv. B%isGlobal)) then
       call dm_printf("Error in dm_yjoin: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    call mat_yjoin(A%x,A%nx,A%ny,A%nz,B%x,B%nx,B%ny,B%nz,C%x,ierr)
    call dm_set_implicit(C,ierr)
    C%isGlobal = A%isGlobal
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
       print *, "Error in dm_yjoin: Matrix A and Matrix B &
            &should have the same distribution."
       print*, A%nx,A%ny,A%nz,B%nx,B%ny,B%nz
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
       print *, "Error in dm_mult: the column of A matrix &
            &should equal to the row of B matrix, and the &
            &number of z dimension should be same."
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
    B%nz=A%nz	
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
  end function dm_mult3

  function dm_mult4(alpha,A) result(B)
    implicit none
    real,           intent(in)	::  alpha 
    type(Matrix),	intent(in)	::  A 
    type(Matrix)                ::	B
    integer						::	ierr
    B=dm_mult2(real(alpha,8),A)
    call dm_set_implicit(B,ierr)
  end function dm_mult4

  function dm_mult5(A,alpha) result(B)
    implicit none
    real,           intent(in)	::  alpha 
    type(Matrix),	intent(in)	::  A 
    type(Matrix)                ::	B
    integer						::	ierr
    B=dm_mult2(real(alpha,8),A)
    call dm_set_implicit(B,ierr)
  end function dm_mult5

  function dm_mult6(A,alpha) result(B)
    implicit none
    integer,        intent(in)	::  alpha 
    type(Matrix),	intent(in)	::  A 
    type(Matrix)                ::	B
    integer						::	ierr
    B=dm_mult2(real(alpha,8),A)
    call dm_set_implicit(B,ierr)
  end function dm_mult6

  function dm_mult7(alpha,A) result(B)
    implicit none
    integer,        intent(in)	::  alpha 
    type(Matrix),	intent(in)	::  A 
    type(Matrix)                ::	B
    integer						::	ierr
    B=dm_mult2(real(alpha,8),A)
    call dm_set_implicit(B,ierr)
  end function dm_mult7

  ! -----------------------------------------------------------------------
  !> C=A.*B element-wise matrix-matrix multiplication
  ! -----------------------------------------------------------------------
  function dm_emult(A,B) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer						::	ierr
    !print *, "(A%nx,A%ny,A%nz)=",A%nx,A%ny,A%nz,"(B%nx,B%ny,B%nz)=",B%nx,B%ny,B%nz 
    if((A%nx/=B%nx) .or. (A%ny/=B%ny) .or. (A%nz/=B%nz) &
         .or.  (A%isGlobal .neqv. B%isGlobal)) then
       print *, "Error in dm_emult: the A matrix &
            &and B matrix should have the same distribution."
       call abort()
       !call backtrace()
       stop	
    endif

    call mat_emult(A%x,A%nx,A%ny,A%nz,B%x,B%nx,B%ny,B%nz,C%x,ierr)
    C%isGlobal=A%isGlobal
    C%nx=A%nx
    C%ny=A%ny
    C%nz=A%nz
    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif
  end function dm_emult


  ! -----------------------------------------------------------------------
  !> C=A./B Element-wise matrix-matrix divide
  ! -----------------------------------------------------------------------
  function dm_ediv(A,B) result(C)
    implicit none
    type(Matrix),intent(in)	::  A,B 
    type(Matrix)              	::	C
    integer			::	ierr

    if((A%nx/=B%nx) .or. (A%ny/=B%ny) .or. (A%nz/=B%nz) .or. &
         (A%isGlobal .neqv. B%isGlobal)) then
       print *, "Error in dm_ediv: the A matrix and B &
            &matrix should have the same distribution." 
       stop	
    endif

    call mat_ediv(A%x,A%nx,A%ny,A%nz,B%x,B%nx,B%ny,B%nz,C%x,ierr)
    C%isGlobal=A%isGlobal
    C%nx=A%nx
    C%ny=A%ny
    C%nz=A%nz
    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif
  end function dm_ediv


  ! -----------------------------------------------------------------------
  ! B=repmat(A,m,n,k)
  ! -----------------------------------------------------------------------
  function dm_rep(A,m,n,k) result(B)
    implicit none
    type(Matrix),   intent(in)	::  A 
    integer,	    intent(in)	::  m,n,k 
    type(Matrix)              	::	B
    integer						::	ierr
    call mat_rep(A%x,A%nx,A%ny,A%nz,m,n,k,B%x,ierr)
    B%isGlobal=A%isGlobal
    B%nx=m*A%nx
    B%ny=n*A%ny
    B%nz=k*A%nz
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function dm_rep


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
    type(Matrix),   intent(in)	::  A 
    integer,	    intent(in)	::  ndim 
    type(Matrix)              	::	B
    integer						::	ierr
    call mat_sum(A%x,A%nx,A%ny,A%nz,ndim,B%x,ierr)
    B%isGlobal=A%isGlobal

    if(ndim == 1) then
       B%nx = 1
       B%ny = A%ny
       B%nz = A%nz
    else if ( ndim == 2 ) then
       B%nx = A%nx
       B%ny = 1
       B%nz = A%nz
    else 
       print*, "Error: ndim can be only 1 or 2"
       stop
    endif

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
    
    call dm_set_implicit(B,ierr)
  end function 

  function dm_sum_all(A) result(val)
    implicit none
    type(Matrix),   intent(in)	:: A 
    real(kind=8)  :: val
    integer			:: ierr

    call mat_sum_all(A%x, A%nx, A%ny, A%nz, val, ierr)
    
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
    Y%nx = X%nx
    Y%ny = X%ny
    Y%nz = X%nz    
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
    Y%nx = X%nx
    Y%ny = X%ny
    Y%nz = X%nz    
    call dm_set_explicit(Y,ierr)

    if (X%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(X%x,ierr)
    endif
  end subroutine 

  subroutine dm_axpy3(Y,a,X,ierr)
    implicit none
    type(Matrix),   intent(in)		::  X 
    integer,        intent(in)		::  a
    type(Matrix), 	intent(inout) 	::  Y 
    integer,		intent(out)	::	ierr
    call mat_axpy(Y%x,real(a,kind=8),X%x,ierr)
    Y%isGlobal=X%isGlobal
    Y%nx = X%nx
    Y%ny = X%ny
    Y%nz = X%nz    
   
    call dm_set_explicit(Y,ierr)

    if (X%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(X%x,ierr)
    endif
  end subroutine 


  ! -----------------------------------------------------------------------
  !> Compute Y = a*Y + X.
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
  !> B = A^T. The transpose only performed on xy-plane.
  !! i.e. (nx,ny,nz)->(ny,nx,nz)
  ! -----------------------------------------------------------------------
  function dm_trans(A) result(B)
    implicit none
    type(Matrix),	intent(in) ::  A 
    type(Matrix)              	::	B
    integer						::	ierr

    call mat_trans(A%x,B%x,ierr)
    B%isGlobal=A%isGlobal
    B%nx = A%ny
    B%ny = A%nx
    B%nz = A%nz
    
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function dm_trans


  ! -----------------------------------------------------------------------
  !> C = A*B^T
  ! -----------------------------------------------------------------------
  function dm_xyt(A,B) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer						::	ierr

    if(A%nx .ne. B%ny .or. A%nz.ne.B%nz) then
       call dm_printf("Error: dimension of A and B &
            &do not match for multiplication.", ierr)
       stop
    endif
    
    if(A%isGlobal .neqv. B%isGlobal) then
       call dm_printf("Error in dm_xyt: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif
    
    call mat_xyt(A%x,B%x,C%x,ierr)
    C%isGlobal=A%isGlobal
    C%nx = A%nx
    C%ny = B%nx
    C%nz = B%nz
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
    C%nx = A%ny
    C%ny = B%ny
    C%nz = A%nz
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
    B%nx = A%nx
    B%ny = A%ny
    B%nz = A%nz
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function dm_exp


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
    B%nx = A%nx
    B%ny = A%ny
    B%nz = A%nz
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function dm_log

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
    B%nx = A%nx
    B%ny = A%ny
    B%nz = A%nz
    call dm_set_implicit(B,ierr)
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function dm_squ


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
    B%nx = A%nx
    B%ny = A%ny
    B%nz = A%nz
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function dm_cube



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
    B%nx = A%nx
    B%ny = A%ny
    B%nz = A%nz
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function dm_sqrt

  ! -----------------------------------------------------------------------
  ! B=abs(A) 
  ! -----------------------------------------------------------------------
  function dm_abs(A) result(B)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A 
    type(Matrix)              	::	B
    integer						::	ierr
    call mat_math(A%x,MAT_MATH_ABS,B%x,ierr)
    B%isGlobal=A%isGlobal
    B%nx = A%nx
    B%ny = A%ny
    B%nz = A%nz
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function dm_abs

  ! -----------------------------------------------------------------------
  ! B=pow(A, p) 
  ! -----------------------------------------------------------------------
  function dm_pow1(A, p) result(B)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A
    real(kind=8), intent(in) :: p
    type(Matrix)              	::	B
    integer						::	ierr
    call mat_math(A%x,MAT_MATH_POW,B%x,ierr, p)
    B%isGlobal=A%isGlobal
    B%nx = A%nx
    B%ny = A%ny
    B%nz = A%nz
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function

  function dm_pow2(A, p) result(B)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A
    real(kind=4), intent(in) :: p
    type(Matrix)        	::	B
    integer	::	ierr

    B = dm_pow1(A, real(p, kind=8))
    
    call dm_set_implicit(B,ierr)
  end function

  function dm_pow3(A, p) result(B)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A
    integer, intent(in) :: p
    type(Matrix)        	::	B
    integer	::	ierr

    B = dm_pow1(A, real(p, kind=8))
    
    call dm_set_implicit(B,ierr)
  end function
  
  ! -----------------------------------------------------------------------
  ! Solve Ax=b 
  ! -----------------------------------------------------------------------
  function dm_solve(A,B) result(X)
    implicit none
    type(Matrix),intent(in)  :: A 
    type(Matrix),intent(in)  ::	B
    type(Matrix)             ::	X
    integer		     ::	ierr

    if(A%isGlobal .neqv. B%isGlobal) then
       call dm_printf("Error in dm_solve: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    if(A%nx.ne.B%nx .or. A%nz.ne.B%nz .or. B%ny.ne.1) then
       call dm_printf("Error in dm_solve: matrix size error for A, B and X",ierr)
       stop
    endif

    call mat_solve(A%x,B%x,X%x,A%nx,A%ny,A%nz,ierr)
    X%isGlobal=A%isGlobal
    X%nx=A%nx
    X%ny=1
    X%nz=A%nz
    call dm_set_implicit(X,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif
  end function dm_solve

  
  ! -----------------------------------------------------------------------
  ! Load a standard row-cloumn file into a matrix 
  ! -----------------------------------------------------------------------
  subroutine dm_load(filename,varname, A, isGlobal, ierr)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname
    logical, intent(in) :: isGlobal 
    type(Matrix), intent(out) :: A 
    integer, intent(out) :: ierr
    logical :: file_exists = .false.
    
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
       print*, "Error: File ",filename," does not exists."
       call dm_abort(ierr)
    endif

    call mat_load(filename, varname, A%x,A%nx,A%ny,A%nz, isGlobal, ierr)
    call dm_set_explicit(A, ierr)
  end subroutine

  subroutine dm_load3d(filename,varname, A, isGlobal, ierr)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname
    logical, intent(in) :: isGlobal 
    type(Matrix), intent(out) :: A 
    integer, intent(out) :: ierr
    logical :: file_exists

    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
       print*, "Error: File ",filename," does not exists."
       stop
    endif
    
    call mat_load3d(filename, varname, A%x,A%nx,A%ny,A%nz, isGlobal, ierr)
    call dm_set_explicit(A, ierr)
  end subroutine


  subroutine dm_save(filename, varname, A, ierr)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname
    type(Matrix) :: A 
    integer, intent(out) :: ierr

    call mat_save(filename, varname, A%x,A%nx,A%ny,A%nz,A%isGlobal,ierr)
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end subroutine 

  subroutine dm_save3d(filename, varname, A, ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname
    type(Matrix) :: A 
    integer, intent(out) :: ierr
    integer :: size
    
    call mat_save3d(filename, varname, A%x,A%nx,A%ny,A%nz,A%isGlobal, ierr)
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end subroutine 
  
  subroutine dm_setdiag1(A, val, ierr)
    implicit none
    type(Matrix), intent(inout) :: A
    real(kind=8) :: val
    integer, intent(inout) :: ierr
    
    call mat_setdiag(A%x, val, A%nx, A%ny, A%nz, ierr)

    if(A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x, ierr)
    endif
  end subroutine 

  subroutine dm_setdiag2(A, val, ierr)
    implicit none
    type(Matrix), intent(inout) :: A
    real(kind=4) :: val
    integer, intent(inout) :: ierr
    
    call dm_setdiag1(A, real(val,kind=8), ierr)
  end subroutine 

  subroutine dm_setdiag3(A, val, ierr)
    implicit none
    type(Matrix), intent(inout) :: A
    integer :: val
    integer, intent(inout) :: ierr
    
    call dm_setdiag1(A, real(val,kind=8), ierr)
  end subroutine 

  function dm_gendiag1(m, n, k, val, isGlobal) result(C)
    implicit none
    type(Matrix) :: C
    logical, intent(in), optional :: isGlobal
    integer :: ierr
    integer, intent(in) :: m, n, k
    real(kind=8), intent(in) :: val

    if(present(isGlobal)) then
       !call dm_zeros(C, m, n, k, isGlobal, ierr)
       C = dm_zeros(m, n, k, isGlobal)       
    else
       !call dm_zeros(C, m, n, k, .true., ierr)
       C = dm_zeros(m, n, k)              
    endif

    call dm_setdiag(C, val, ierr)
    call dm_set_implicit(C, ierr)    
  end function 

  function dm_gendiag2(m, n, k, val, isGlobal) result(C)
    implicit none
    type(Matrix) :: C
    logical, intent(in), optional :: isGlobal
    integer :: ierr
    integer, intent(in) :: m, n, k    
    real, intent(in) :: val

    C = dm_gendiag1(m, n, k, real(val, kind=8), isGlobal)
    call dm_set_implicit(C, ierr)        
  end function
  
  function dm_gendiag3(m, n, k, val, isGlobal) result(C)
    implicit none
    type(Matrix) :: C
    logical, intent(in), optional :: isGlobal
    integer :: ierr
    integer, intent(in) :: m, n, k    
    integer, intent(in) :: val

    C = dm_gendiag1(m, n, k, real(val, kind=8), isGlobal)    
    call dm_set_implicit(C, ierr)
  end function

  function dm_getdiag(A) result(C)
    implicit none
    type(Matrix), intent(in) :: A
    type(Matrix) :: C
    integer :: ierr
    
    call mat_getdiag(A%x, C%x, A%nx, A%ny, A%nz, ierr)
    C%nx = A%nx
    C%ny = 1
    C%nz = A%nz
    C%isGlobal = A%isGlobal
    
    if(A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x, ierr)
    endif
  end function 
  
  ! -----------------------------------------------------------------------
  ! A(m,n)=value 
  ! -----------------------------------------------------------------------
  subroutine dm_setvalue1(A,m,n,k,value,ierr)
    implicit none
    type(Matrix)	            ::  A 
    integer,	    intent(in)	::	m,n,k
    real(kind=8),   intent(in)	::	value
    integer,		intent(out)	::	ierr

    call mat_setvalue(A%x,A%nx,A%ny,A%nz,m,n,k,value,ierr)
  end subroutine 

  subroutine dm_setvalue2(A,m,n,k,value,ierr)
    implicit none
    type(Matrix)	            ::  A 
    integer,	    intent(in)	::	m,n,k
    integer,    	intent(in)	::	value
    integer,		intent(out)	::	ierr

    call mat_setvalue(A%x,A%nx,A%ny,A%nz,m,n,k,real(value,8),ierr)
  end subroutine 

  subroutine dm_setvalue3(A,m,n,k,value,ierr)
    implicit none
    type(Matrix)	            ::  A 
    integer,	    intent(in)	::	m,n,k
    real,    		intent(in)	::	value
    integer,		intent(out)	::	ierr

    call mat_setvalue(A%x,A%nx,A%ny,A%nz,m,n,k,real(value,8),ierr)
  end subroutine 


  ! -----------------------------------------------------------------------
  ! B=A(rows,cols). Get sub matrix.
  ! -----------------------------------------------------------------------
  function dm_getsub(A,idxm,idxn,idxk) result(B)
    implicit none
    type(Matrix),	intent(in)	::	A
    integer,		intent(in)	::	idxm(:),idxn(:),idxk(:)
    type(Matrix)			::	B
    integer				::	ierr

    !print *, "A%isGLobal=",A%isGlobal
    call mat_getsub(A%x,A%nx,A%ny,A%nz,idxm,idxn,idxk,B%x,ierr)
    B%isGLobal=A%isGlobal
    B%nx = size(idxm)
    B%ny = size(idxn)
    B%nz = size(idxk)
    
    call dm_set_implicit(B,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end function 


  ! -----------------------------------------------------------------------
  ! B=A(:,i). Get a column from A.
  ! -----------------------------------------------------------------------
  ! function dm_getcol(A,n) result(B)
  !   implicit none
  !   type(Matrix),	intent(in)	::	A
  !   integer,        intent(in)  ::  n
  !   type(Matrix)				::	B
  !   integer						::	ierr
  !   integer 					:: 	i

  !   call mat_getsub(A%x, (/(i,i=0,A%nrow-1)/), (/n/), B%x, ierr)
  !   B%isGlobal=A%isGlobal
  !   call dm_set_implicit(B,ierr)

  !   if (A%xtype==MAT_XTYPE_IMPLICIT) then
  !      call mat_destroy(A%x,ierr)
  !   endif
  ! end function dm_getcol


  ! -----------------------------------------------------------------------
  ! B=A(m,:). Get a row from A.
  ! -----------------------------------------------------------------------
  ! function dm_getrow(A,n) result(B)
  !   implicit none
  !   type(Matrix),	intent(in)	::	A
  !   integer,        intent(in)  ::  n
  !   type(Matrix)				::	B
  !   integer						::	ierr
  !   integer,allocatable			::	rows(:),cols(:)
  !   integer 					:: 	i

  !   allocate(rows(1),cols(A%ncol))	
  !   rows(1)=n
  !   do i=1,A%ncol
  !      cols(i)=i-1
  !   enddo

  !   call mat_getsub(A%x, A%nx,A%ny,A%nz,rows, cols, B%x, ierr)
  !   B%isGlobal=A%isGlobal
  !   call dm_set_implicit(B,ierr)

  !   deallocate(rows,cols)	

  !   if (A%xtype==MAT_XTYPE_IMPLICIT) then
  !      call mat_destroy(A%x,ierr)
  !   endif
  ! end function dm_getrow


  ! -----------------------------------------------------------------------
  ! Set local values in A.
  ! -----------------------------------------------------------------------
  subroutine dm_setvalues1(A,idxm,idxn,idxk,v,ierr)
    implicit none
    type(Matrix),	intent(in)	::	A
    integer,		intent(in)	::	idxm(:),idxn(:),idxk(:)
    real(kind=8),	intent(in)	::	v(:)	
    integer,		intent(out)	::	ierr

    call mat_setvalues(A%x,A%nx,A%ny,A%nz,idxm,idxn,idxk,v,ierr) 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end subroutine dm_setvalues1

  subroutine dm_setvalues2(A,idxm,idxn,idxk,v,ierr)
    implicit none
    type(Matrix),	intent(in)	::	A
    integer,		intent(in)	::	idxm(:),idxn(:),idxk(:)
    real,		intent(in)	::v(:)	
    integer,		intent(out)	::	ierr

    call mat_setvalues(A%x,A%nx,A%ny,A%nz,idxm,idxn,idxk,real(v,8),ierr) 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end subroutine dm_setvalues2

  subroutine dm_setvalues3(A,idxm,idxn,idxk,v,ierr)
    implicit none
    type(Matrix),	intent(in)	::	A
    integer,		intent(in)	::	idxm(:),idxn(:),idxk(:)
    integer,		intent(in)	::	v(:)	
    integer,		intent(out)	::	ierr

    call mat_setvalues(A%x,A%nx,A%ny,A%nz,idxm,idxn,idxk,real(v,8),ierr) 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end subroutine dm_setvalues3


  ! -----------------------------------------------------------------------
  ! Get local values in A.
  ! -----------------------------------------------------------------------
  subroutine dm_getvalues(A,idxm,idxn,idxk,v,ierr)
    implicit none
    type(Matrix),	intent(in)	::	A
    integer,		intent(in)	::	idxm(:),idxn(:),idxk(:)
    real(kind=8),	intent(inout)	::	v(:)	
    integer,		intent(out)	::	ierr

    call mat_getvalues(A%x,A%isGlobal,A%nx,A%ny,A%nz,idxm,idxn,idxk,v,ierr) 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif
  end subroutine dm_getvalues

  ! -----------------------------------------------------------------------
  ! Get a single value in A.
  ! -----------------------------------------------------------------------
  function dm_getvalue(A, idxm, idxn, idxk) result(res)
    implicit none
    type(Matrix),	intent(in)	::	A
    integer,		intent(in)	::	idxm,idxn,idxk
    real(kind=8) :: res
    real(kind=8)		::	v(1)	
    integer		        ::	ierr

    call dm_getvalues(A, (/idxm/),(/idxn/),(/idxk/),v,ierr) 
    res = v(1)
  end function
  
  ! ! -----------------------------------------------------------------------
  ! ! Get local values in A.
  ! ! -----------------------------------------------------------------------
  ! subroutine dm_getvalues2(A,idxm,idxn,idxk,v,ierr)
  !   implicit none
  !   type(Matrix),	intent(in)	::	A
  !   integer,		intent(in)	::	idxm(:),idxn(:),idxk(:)
  !   real(kind=8),	intent(inout)	::	v(:)	
  !   integer,		intent(out)	::	ierr

  !   call mat_getvalues2(A%x,A%isGlobal,A%nx,A%ny,A%nz,idxm,idxn,idxk,v,ierr) 

  !   if (A%xtype==MAT_XTYPE_IMPLICIT) then
  !      call mat_destroy(A%x,ierr)
  !   endif
  ! end subroutine dm_getvalues


  ! -----------------------------------------------------------------------
  !> Norm(A)
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
    type(Matrix),intent(in) :: A
    real(kind=8)            :: res	
    integer                 :: ierr
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

    if(A%nx.ne.B%nx .or. A%ny.ne.B%ny .or. A%nz.ne.A%nz) then
       print *, "Error: Matrix A and matrix B should have the same size"
       stop	
    endif

    if(A%isGlobal .neqv. B%isGlobal) then
       call dm_printf("Error in dm_lt: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_LT,C%x,A%nx,A%ny,A%nz,ierr)

    C%nx = A%nx
    C%ny = A%ny
    C%nz = A%nz
    C%isGlobal=A%isGlobal
    
    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif

    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif

  end function dm_lt1

  function dm_lt2(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	B
    type(Matrix)                ::	C
    integer::	ierr

    B=dm_constants(A%nx,A%ny,A%nz,alpha, A%isGlobal)
    C=dm_lt1(A,B)
    call dm_destroy(B, ierr)
    call dm_set_implicit(C,ierr)
  end function dm_lt2

  function dm_lt3(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,	        intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_lt2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_lt3

  function dm_lt4(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,	    intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr
    
    C=dm_lt2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_lt4

  ! -----------------------------------------------------------------------
  !> C= A<=B
  ! -----------------------------------------------------------------------
  function dm_le1(A,B) result(C)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer						::	ierr

    if(A%nx.ne.B%nx .or. A%ny.ne.B%ny .or. A%nz.ne.A%nz) then
       print *, "Error: Matrix A and matrix B should have the same size"
       stop	
    endif

    if(A%isGlobal .neqv. B%isGlobal) then
       call dm_printf("Error in dm_le: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_LE,C%x,A%nx,A%ny,A%nz,ierr)

    C%nx = A%nx
    C%ny = A%ny
    C%nz = A%nz
    C%isGlobal=A%isGlobal

    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif

    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif

  end function dm_le1

  function dm_le2(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	B
    type(Matrix)                ::	C
    integer::	ierr

    B=dm_constants(A%nx,A%ny,A%nz,alpha, A%isGlobal)
    C=dm_le1(A,B)
    call dm_destroy(B, ierr)
    call dm_set_implicit(C,ierr)
  end function dm_le2

  function dm_le3(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,	        intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_le2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_le3

  function dm_le4(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,	    intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_le2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_le4

  ! -----------------------------------------------------------------------
  !> C= A>B
  ! -----------------------------------------------------------------------
  function dm_gt1(A,B) result(C)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer			::	ierr

    if(A%nx.ne.B%nx .or. A%ny.ne.B%ny .or. A%nz.ne.A%nz) then
       print *, "Error: Matrix A and matrix B should have the same size"
       stop	
    endif

    if(A%isGlobal .neqv. B%isGlobal) then
       call dm_printf("Error in dm_gt: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_GT,C%x,A%nx,A%ny,A%nz,ierr)

    C%nx = A%nx
    C%ny = A%ny
    C%nz = A%nz
    C%isGlobal=A%isGlobal

    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif

    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif

  end function dm_gt1

  function dm_gt2(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	B
    type(Matrix)                ::	C
    integer::	ierr

    B=dm_constants(A%nx,A%ny,A%nz,alpha, A%isGlobal)
    C=dm_gt1(A,B)
    call dm_destroy(B, ierr)
    call dm_set_implicit(C,ierr)
  end function dm_gt2

  function dm_gt3(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,	        intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_gt2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_gt3

  function dm_gt4(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,	    intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_gt2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_gt4

  ! -----------------------------------------------------------------------
  !> C= A>=B
  ! -----------------------------------------------------------------------
  function dm_ge1(A,B) result(C)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer			::	ierr

    if(A%nx.ne.B%nx .or. A%ny.ne.B%ny .or. A%nz.ne.A%nz) then
       print *, "Error: Matrix A and matrix B should have the same size"
       stop	
    endif

    if(A%isGlobal .neqv. B%isGlobal) then
       call dm_printf("Error in dm_ge: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_GE,C%x,A%nx,A%ny,A%nz,ierr)

    C%nx = A%nx
    C%ny = A%ny
    C%nz = A%nz
    C%isGlobal=A%isGlobal

    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif

    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif

  end function dm_ge1

  function dm_ge2(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	B
    type(Matrix)                ::	C
    integer::	ierr

    B=dm_constants(A%nx,A%ny,A%nz,alpha, A%isGlobal)
    C=dm_ge1(A,B)
    call dm_destroy(B, ierr)
    call dm_set_implicit(C,ierr)
  end function dm_ge2

  function dm_ge3(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,	        intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_ge2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_ge3

  function dm_ge4(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,	    intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_ge2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_ge4

  ! -----------------------------------------------------------------------
  !> C= A==B
  ! -----------------------------------------------------------------------
  function dm_eq1(A,B) result(C)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer			::	ierr

    if(A%nx.ne.B%nx .or. A%ny.ne.B%ny .or. A%nz.ne.A%nz) then
       print *, "Error: Matrix A and matrix B should have the same size"
       stop	
    endif

    if(A%isGlobal .neqv. B%isGlobal) then
       call dm_printf("Error in dm_eq: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_EQ,C%x,A%nx,A%ny,A%nz,ierr)

    C%nx = A%nx
    C%ny = A%ny
    C%nz = A%nz
    C%isGlobal=A%isGlobal

    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif

    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif

  end function dm_eq1

  function dm_eq2(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	B
    type(Matrix)                ::	C
    integer::	ierr

    B=dm_constants(A%nx,A%ny,A%nz,alpha, A%isGlobal)
    C=dm_eq1(A,B)
    call dm_destroy(B, ierr)
    call dm_set_implicit(C,ierr)
  end function dm_eq2

  function dm_eq3(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,	        intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_eq2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_eq3

  function dm_eq4(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,	    intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_eq2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_eq4

  ! -----------------------------------------------------------------------
  !> C= (A/=B)
  ! -----------------------------------------------------------------------
  function dm_nq1(A,B) result(C)
    implicit none
#include "mat_type.h"
    type(Matrix),	intent(in)	::  A 
    type(Matrix),	intent(in)	::  B 
    type(Matrix)              	::	C
    integer			::	ierr

    if(A%nx.ne.B%nx .or. A%ny.ne.B%ny .or. A%nz.ne.A%nz) then
       print *, "Error: Matrix A and matrix B should have the same size"
       stop	
    endif

    if(A%isGlobal .neqv. B%isGlobal) then
       call dm_printf("Error in dm_nq: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    call mat_compare(A%x,B%x,MAT_COMPARE_NQ,C%x,A%nx,A%ny,A%nz,ierr)

    C%nx = A%nx
    C%ny = A%ny
    C%nz = A%nz
    C%isGlobal=A%isGlobal

    call dm_set_implicit(C,ierr)

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x,ierr)
    endif

    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif

  end function dm_nq1

  function dm_nq2(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real(kind=8),	intent(in)	::  alpha 
    type(Matrix)                ::	B
    type(Matrix)                ::	C
    integer::	ierr

    B=dm_constants(A%nx,A%ny,A%nz,alpha, A%isGlobal)
    C=dm_nq1(A,B)
    call dm_destroy(B, ierr)
    call dm_set_implicit(C,ierr)
  end function dm_nq2

  function dm_nq3(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    real,	        intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_nq2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_nq3

  function dm_nq4(A,alpha) result(C)
    implicit none
    type(Matrix),	intent(in)	::  A 
    integer,	    intent(in)	::  alpha 
    type(Matrix)                ::	C
    integer::	ierr

    C=dm_nq2(A,real(alpha, kind=8))
    call dm_set_implicit(C,ierr)
  end function dm_nq4


  ! -----------------------------------------------------------------------
  ! sparse(i,j,A,B,ierr)
  ! -----------------------------------------------------------------------
  function dm_sparse(Ind_m,Ind_n,Ind_k,A,m,n,k) result(B)
    implicit none
    type(Matrix), intent(in)	::  Ind_m,Ind_n,Ind_k 
    type(Matrix), intent(in)	::  A 
    integer,	  intent(in)	::  m,n,k
    logical :: isGlobal
    
    type(Matrix)		::  B
    integer			::  ierr

    if((Ind_m%isGlobal .neqv. Ind_n%isGlobal) &
         .or. (A%isGlobal .neqv. Ind_m%isGlobal) &
         .or. (A%isGlobal .neqv. Ind_k%isGlobal)) then
       call dm_printf("Error in dm_sparse: Matrix Ind_m, Ind_n, Ind_k and A &
            &should have the same distribution.",ierr)
       stop
    endif

    if(Ind_m%nz.ne.1 .or. Ind_n%nz.ne.1 .or. Ind_k%nz.ne.1 &
         .or. Ind_m%ny.ne.1 .or. Ind_n%ny.ne.1 .or. Ind_k%ny.ne.1 &
         .or. Ind_m%nx.ne.Ind_n%nx .or. Ind_m%nx.ne.Ind_k%nx ) then
       call dm_printf("Error (dm_sparse): Ind_m, Ind_n and Ind_k &
            &must be a Nx1x1 matrix", ierr)
       stop
    endif

    isGlobal = A%isGlobal
    
    B=dm_zeros(m,n,k, isGlobal)
    call mat_sparse(Ind_m%x,Ind_n%x,Ind_k%x,A%x,m,n,k,B%x,ierr)
    call dm_set_implicit(B,ierr)
  end function dm_sparse



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
  end subroutine dm_cart2sph


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
  end subroutine dm_set_implicit


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
  end subroutine dm_set_explicit


  ! -----------------------------------------------------------------------
  ! Set the diagonal of A to constant value 
  ! -----------------------------------------------------------------------
  ! subroutine dm_setdiag1(A,value,ierr) 
  !   implicit none
  !   type(Matrix),	intent(inout)	::  A 
  !   real(kind=8),	intent(in)		::  value
  !   integer							::	ierr

  !   call mat_setdiag(A%x,value,ierr)
  !   call dm_set_explicit(A,ierr)
  ! end subroutine dm_setdiag1

  ! subroutine dm_setdiag2(A,value,ierr) 
  !   implicit none
  !   type(Matrix),	intent(inout)	::  A 
  !   real,			intent(in)		::  value
  !   integer							::	ierr

  !   call mat_setdiag(A%x,real(value,kind=8),ierr)
  !   call dm_set_explicit(A,ierr)
  ! end subroutine dm_setdiag2

  ! subroutine dm_setdiag3(A,value,ierr) 
  !   implicit none
  !   type(Matrix),	intent(inout)	::  A 
  !   integer,		intent(in)		::  value
  !   integer							::	ierr

  !   call mat_setdiag(A%x,real(value,kind=8),ierr)
  !   call dm_set_explicit(A,ierr)
  ! end subroutine dm_setdiag3


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
       call dm_printf("Error in dm_setcol: Matrix A and B &
            &should have the same distribution.",ierr)
       stop
    endif

    if(A%nrow/= B%nrow)then
       print *, "Error in dm_setcol: Matrix A and Matrix B &
            &should have the same row size."
       stop	
    endif

    call mat_setcol(A%x,idxn,B%x,ierr)	
    call dm_set_explicit(A,ierr)
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x,ierr)
    endif
  end subroutine dm_setcol

  subroutine dm_max(A, val, pos, ierr)
    implicit none
    type(Matrix), intent(in) :: A
    integer, intent(out) :: pos(3)
    integer, intent(out) :: ierr
    real(kind=8), intent(out) :: val
    
    call mat_max_min(A%x, A%nx, A%ny, A%nz, val, &
         pos, .false., ierr)
    
    if(A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x, ierr)
    endif
  end subroutine 
  
  subroutine dm_min(A, val, pos, ierr)
    implicit none    
    type(Matrix), intent(in) :: A
    integer, intent(out) :: pos(3)
    integer, intent(out) :: ierr
    real(kind=8), intent(out) :: val
    
    call mat_max_min(A%x, A%nx, A%ny, A%nz, val, &
         pos, .true., ierr)
    
    if(A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x, ierr)
    endif
  end subroutine

  function dm_trid(A, B, C) result(D)
    implicit none
    type(Matrix), intent(in)  :: A,B,C
    type(Matrix) :: D
    integer :: ierr

    D%isGlobal = A%isGlobal
    D%nx = A%nz
    D%ny = A%nz
    D%nz = A%nx * A%ny
    call mat_trid(A%x, B%x, C%x, A%nx, A%ny, A%nz, D%x, ierr)

    if(A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x, ierr)
    endif
    if(B%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(B%x, ierr)
    endif
    if(C%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(C%x, ierr)
    endif
    call dm_set_implicit(D, ierr)
  end function

  function dm_trid1(A) result(D)
    implicit none
    type(Matrix), intent(in)  :: A
    type(Matrix) :: D
    integer :: ierr

    D%isGlobal = A%isGlobal
    D%nx = A%nz 
    D%ny = 1
    D%nz = A%nx * A%ny
    call mat_trid1(A%x, A%nx, A%ny, A%nz, D%x, ierr)

    if(A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x, ierr)
    endif
    call dm_set_implicit(D, ierr)
  end function

  function dm_trid2(A, Dnx, Dny, Dnz) result(D)
    implicit none
    type(Matrix), intent(in)  :: A
    type(Matrix) :: D
    integer,intent(in) :: Dnx, Dny, Dnz
    integer :: ierr

    D%isGlobal = A%isGlobal
    D%nx = Dnx 
    D%ny = Dny
    D%nz = Dnz
    call mat_trid2(A%x, A%nx, A%ny, A%nz, D%x,Dnx, Dny, Dnz, ierr)

    if(A%xtype==MAT_XTYPE_IMPLICIT) then
       call mat_destroy(A%x, ierr)
    endif
    call dm_set_implicit(D, ierr)
  end function

  function dm_find1(A) result(C)
    implicit none
    type(Matrix), intent(in) :: A
    type(Matrix) :: C
    integer :: ierr
    integer :: nonzero_count

    call mat_find1(A%x, A%nx, A%ny, A%nz, C%x, nonzero_count, ierr)

    if(nonzero_count == 0) then
       C%nx = 0
       C%ny = 0
       C%nz = 0
    else
       C%nx = nonzero_count
       C%ny = 1
       C%nz = 1
    endif
    C%isGlobal = A%isGlobal

    !print*, "nonzero_count=", nonzero_count
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif
    call dm_set_implicit(C, ierr)
  end function dm_find1

  function dm_find2(A) result(C)
    implicit none
    type(Matrix), intent(in) :: A
    integer, allocatable :: C(:)
    integer :: ierr

    call mat_find2(A%x, A%nx, A%ny, A%nz, C, ierr)

    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif
  end function dm_find2

  function dm_inverse(A) result(X)
    implicit none
    type(Matrix), intent(in) :: A
    type(Matrix) :: X
    integer :: ierr
    character(len=200) :: msg
    
    if(A%nx /= A%ny .or. A%nz /= 1) then
       write(msg, *) "dm_inverse: matrix dimension error! size(A)=", &
            A%nx,A%ny,A%nz
       
       call dm_printf(msg, ierr)
       call abort()
    endif
    
    X%x = mat_inverse(A%x, A%nx)
    X%nx = A%nx
    X%ny = A%ny
    X%nz = 1
    call dm_set_implicit(X, ierr)
    
  end function dm_inverse
end module dm
