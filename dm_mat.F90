! -----------------------------------------------------------------------
!> Distributed Matrix Computing Library
! -----------------------------------------------------------------------
module dm_mat 
  implicit none

contains

  ! -----------------------------------------------------------------------
  !> Generate a logical 3D matrix with 2D Mat 
  ! -----------------------------------------------------------------------
  ! For example, suppose A is a 3*2*2 matrix, we will store A as 
  ! A=[x11 x12   0   0 ]
  ![x21 x22   0   0 ]
  ![x31 x32   0   0 ]
  ![  0   0 Y11 Y12 ]
  ![  0   0 Y21 Y22 ]
  ![  0   0 Y31 Y32 ]
  subroutine mat_create(A,m,n,k,isGlobal,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,  intent(in)  ::  m,n,k    
    PetscBool, intent(in)  ::  isGlobal
    Mat,       intent(out) ::  A
    PetscErrorCode,  intent(out) :: ierr
    PetscInt  ::  i               
    PetscLogEvent               ::  ievent

    call PetscLogEventRegister("mat_create",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    ! generate matrix A with size m*n
    if(isGlobal) then
       call MatCreate(PETSC_COMM_WORLD,A,ierr)
    else
       call MatCreate(PETSC_COMM_SELF,A,ierr);
    endif
    
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*k,n*k,ierr)
    call MatSetFromOptions(A,ierr)
    call MatSetUp(A,ierr)
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_create

  subroutine mat_destroy(A,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,   intent(in)::A
    PetscErrorCode,intent(out)::ierr
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_destroy",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    ! destroy matrix A
    call MatDestroy(A,ierr)
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_destroy

  subroutine mat_view(A,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::A
    PetscErrorCode,intent(out)::ierr
    PetscBool::  isGlobal
    ! view matrix A
    call mat_assemble(A,ierr)
    call mat_gettype(A,isGlobal,ierr)
    if(isGlobal) then
       call MatView(A,PETSC_VIEWER_STDOUT_WORLD, ierr)
    else
       call MatView(A,PETSC_VIEWER_STDOUT_SELF, ierr)
    endif
  end subroutine mat_view

  subroutine vec_view(A,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Vec,intent(in)::A
    PetscErrorCode,intent(out)::ierr
    PetscBool::  isGlobal
    ! view matrix A
    call vec_assemble(A,ierr)
    !call mat_gettype(A,isGlobal,ierr)
    if(isGlobal) then
       call VecView(A,PETSC_VIEWER_STDOUT_WORLD, ierr)
    else
       call VecView(A,PETSC_VIEWER_STDOUT_SELF, ierr)
    endif
  end subroutine vec_view

  subroutine mat_assemble(A,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::A
    PetscErrorCode,intent(out)::ierr
    PetscBool::  assembled
    !call MatAssembled(A,assembled,ierr)
    !if(.not. assembled) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    !endif
  end subroutine mat_assemble

    subroutine vec_assemble(A,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Vec,intent(in)::A
    PetscErrorCode,intent(out)::ierr
    PetscBool::  assembled
    !call MatAssembled(A,assembled,ierr)
    !if(.not. assembled) then
    call VecAssemblyBegin(A,ierr)
    call VecAssemblyEnd(A,ierr)
    !endif
  end subroutine vec_assemble


  ! -----------------------------------------------------------------------
  ! A=0 
  ! -----------------------------------------------------------------------
  subroutine mat_zeros(A,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(out)::A
    PetscErrorCode,intent(out)::ierr
    PetscInt:: ista,iend
    !PetscLogEvent            ::  ievent
    !   call PetscLogEventRegister("mat_zeros",0, ievent, ierr)
    !   call PetscLogEventBegin(ievent,ierr)

    call MatZeroEntries(A,ierr)
    !   call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_zeros

  ! -----------------------------------------------------------------------
  ! A=1. 
  ! -----------------------------------------------------------------------
  subroutine mat_constants(A,m,n,k,alpha,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(out)::A
    PetscErrorCode,intent(out)::ierr
    PetscInt,intent(in)::  m,n,k 
    PetscScalar,intent(in)::  alpha
    PetscInt::  ista,iend,ilocal
    PetscInt,allocatable::idxm(:),idxn(:)
    PetscScalar,allocatable::row(:)
    integer :: i,j
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_constants",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)
    ilocal=iend-ista
    allocate(idxm(1),idxn(n),row(n))
    row=alpha
    do i=ista,iend-1
       idxm(1)=i
       do j=0,n-1
          idxn(j+1)=(i/m)*n+j
       enddo
       !print *,"i=",i,"idxm=",idxm,"idxn=",idxn
       call MatSetValues(A,1,idxm,n,idxn,row,INSERT_VALUES,ierr)
    enddo
    deallocate(idxm,idxn,row)
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_constants


  ! -----------------------------------------------------------------------
  !> The eye function is used to generate the simple and complex identity matrixs. 
  !! For example, if A is a 2*4 matrix, we can use mat_eye(A,ierr) to obtain 
  !> A= [1 0 0 0]
  !> [0 1 0 0]
  !> if A is a 4*2 matrix, then mat_eye(A,ierr) will generate
  !> A= [1 0]
  !> [0 1]
  !> [0 0]
  !> [0 0]
  ! -----------------------------------------------------------------------
  subroutine mat_eye(A,m,n,k,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(out)::A
    PetscErrorCode,intent(out)::ierr
    PetscInt,intent(in)::m,n,k
    PetscInt::nmin
    PetscInt::  ista,iend,xpos,ypos
    PetscScalar::row
    integer :: i,j,counter
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_eye",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call MatZeroEntries(A,ierr)
    nmin=min(m,n)
    row=1.0
    call MatGetOwnershipRange(A,ista,iend,ierr)
    xpos=-1
    ypos=-1
    do i=ista,iend-1
       xpos=mod(i,m)
       if(xpos<n) then
          ypos=i/m*n+xpos
          call MatSetValue(A,i,ypos,row,INSERT_VALUES,ierr)
       endif
       !print *,"i=",i,"ista:iend=",ista,iend,"xpos=",xpos,"ypos=",ypos
       xpos=-1
       ypos=-1
    enddo

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_eye

  !---------------------------------------------------------------------
  !>Create matrix A and then assign a sequence of numbers
  !>Output: A = [1, 2, 3]
  !>            [4, 5, 6]
  !---------------------------------------------------------------------
  subroutine mat_seqs(A,m,n,k,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(out)::A
    PetscErrorCode,intent(out)::ierr
    PetscInt,intent(in)::m,n,k
    PetscInt::nmin
    PetscInt::  ista,iend
    PetscInt, allocatable :: idxn(:), idxm(:)
    PetscScalar, allocatable ::row(:)
    integer :: ik
    PetscLogEvent            ::  ievent
    PetscInt :: ncol, nrow
    integer :: i, j
    
    call PetscLogEventRegister("mat_seqs",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)

    allocate(idxn(n), row(n), idxm(1))

    do i=ista,iend-1
       ik = i / m
       do j=1,n
          idxn(j) = j + ik * n - 1
          row(j) = j+mod(i,m)*n + ik*m*n - 1
       end do
       idxm(1) = i
       call MatSetValues(A, 1, idxm, n, idxn, row, INSERT_VALUES, ierr)
    enddo
    
    deallocate(idxn, row, idxm)
    
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_seqs

  !---------------------------------------------------------------------
  !> Create a random matrix A 
  !---------------------------------------------------------------------
  subroutine mat_rand(A,m,n,k,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(out)::A
    PetscErrorCode,intent(out)::ierr
    PetscInt,intent(in)::m,n,k
    PetscInt::nmin
    PetscInt::  ista,iend
    PetscInt, allocatable :: idxn(:), idxm(:)
    PetscScalar, allocatable ::row(:)
    integer :: ik
    PetscLogEvent            ::  ievent
    PetscInt :: ncol, nrow
    integer :: i, j
    
    call PetscLogEventRegister("mat_seqs",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)

    allocate(idxn(n), row(n), idxm(1))

    do i=ista,iend-1
       ik = i / m
       do j=1,n
          idxn(j) = j + ik * n - 1
          call RANDOM_NUMBER(row)
       end do
       idxm(1) = i
       call MatSetValues(A, 1, idxm, n, idxn, row, INSERT_VALUES, ierr)
    enddo
    
    deallocate(idxn, row, idxm)
    
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_rand

  ! -----------------------------------------------------------------------
  ! The vertical eye plus zero matrix
  ! A= [1 0 0]
  ! [0 1 0]
  ! [0 0 1]
  ! [0 0 0]
  ! [0 0 0]
  ! -----------------------------------------------------------------------
  subroutine mat_veyezero(A,nrow1,nrow2,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(out)::A
    PetscErrorCode,intent(out)::ierr
    PetscInt,       intent(in)::nrow1,nrow2
    PetscInt::  ista,iend
    PetscScalar::row
    integer :: i,j
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_veyezero",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_zeros(A,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)

    do i=ista,iend-1
       do j=0,nrow1-1
          if(i==j) then
             row=real(1.0,kind=8)
             call MatSetValues(A,1,i,1,j,row,INSERT_VALUES,ierr)
          endif
       enddo
    enddo
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_veyezero


  ! -----------------------------------------------------------------------
  ! The vertical zero plus eye matrix
  ! A= [0 0 0]
  ! [0 0 0]
  ! [1 0 0]
  ! [0 1 0]
  ! [0 0 1]
  ! -----------------------------------------------------------------------
  subroutine mat_vzeroeye(A,nrow1,nrow2,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(out)::A
    PetscErrorCode,intent(out)::ierr
    PetscInt,       intent(in)::nrow1,nrow2
    PetscInt::  ista,iend
    PetscScalar::row
    integer :: i,j
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_vzeroeye",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_zeros(A,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)

    do i=ista,iend-1
       do j=0,nrow2-1
          if(i==j+nrow1) then
             row=real(1.0,kind=8)
             call MatSetValues(A,1,i,1,j,row,INSERT_VALUES,ierr)
          endif
       enddo
    enddo
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_vzeroeye



  ! -----------------------------------------------------------------------
  ! B=A 
  ! -----------------------------------------------------------------------
  subroutine mat_copy(A,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::A
    Mat,intent(out)::B
    PetscErrorCode,intent(out)::ierr
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_copy",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    call mat_assemble(A,ierr)
    call MatDuplicate(A,MAT_COPY_VALUES,B,ierr)
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_copy


  ! -----------------------------------------------------------------------
  ! C=[A B] 
  ! -----------------------------------------------------------------------
  subroutine mat_yjoin(A,m1,n1,k1,B,m2,n2,k2,C,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,        intent(in)::  A,B 
    PetscInt,    intent(in)  ::  m1,n1,k1,m2,n2,k2
    Mat,            intent(out)::C
    PetscErrorCode,intent(out)::ierr
    PetscInt::nrow1,ncol1,nrow2,ncol2
    PetscInt::col1,col2,m,n
    PetscInt,allocatable::idxn1(:),idxn2(:),idxn3(:)
    PetscScalar,allocatable::row1(:),row2(:),row3(:)
    PetscInt::  ista,iend
    PetscBool::isGlobal
    integer::i
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_xjoin",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call MatGetSize(A,nrow1,ncol1,ierr)
    call MatGetSize(B,nrow2,ncol2,ierr)
    if(nrow1/=nrow2)then
       print *, "Error in mat_xjoin: Matrix A and Matrix B should have the same row size"
       stop
    endif

    call mat_assemble(A,ierr)
    call mat_assemble(B,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)
    call mat_gettype(A,isGlobal,ierr)
    call mat_create(C,m1,n1+n2,k1,isGlobal,ierr)
    call mat_zeros(C,ierr)
    do i=ista,iend-1
       call MatGetRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
       m=col1
       call MatRestoreRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)

       call MatGetRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
       n=col2
       call MatRestoreRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)

       allocate(idxn1(m),row1(m))
       allocate(idxn2(n),row2(n))
       allocate(idxn3(m+n),row3(m+n))

       call MatGetRow(A,i,col1,idxn1,row1,ierr)
       idxn3(1:m)=i/m1*(n1+n2)+mod(idxn1,n1)
       row3(1:m)=row1
       call MatRestoreRow(A,i,col1,idxn1,row1,ierr)

       call MatGetRow(B,i,col2,idxn2,row2,ierr)
       idxn3((m+1):(m+n))=i/m1*(n1+n2)+n1+mod(idxn2,n2)
       row3((m+1):(m+n))=row2
       call MatRestoreRow(B,i,col2,idxn2,row2,ierr)

       !print *,">i=",i,"idxn1=",idxn1,"idxn2=",idxn2,"idxn3=",idxn3,"row3=",row3
       call MatSetValues(C,1,i,(m+n),idxn3,row3,INSERT_VALUES,ierr)
       deallocate(idxn1,idxn2,idxn3,row1,row2,row3)
    enddo
    call PetscLogEventEnd(ievent,ierr)
  end subroutine 

  ! -----------------------------------------------------------------------
  ! C=[A] 
  !   [B] 
  ! -----------------------------------------------------------------------

  subroutine mat_xjoin(A,m1,n1,k1,B,m2,n2,k2,C,ierr)
    implicit none 
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,            intent(in)  ::  A,B  
    PetscInt,    intent(in)::  m1,n1,k1,m2,n2,k2 
    Mat,            intent(out) ::  C
    PetscErrorCode, intent(out) ::  ierr 
    Mat                         ::  W1,W2,W3
    PetscBool    ::  isGlobal 
    PetscLogEvent               ::  ievent
    call PetscLogEventRegister("mat_yjoin",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    !call mat_assemble(A,ierr)
    !call mat_assemble(B,ierr)

    call mat_trans(A,W1,ierr)
    call mat_trans(B,W2,ierr)
    !print *, "W1="
    !call mat_view(W1,ierr)
    !print *, "W2="
    !call mat_view(W2,ierr)
    call mat_gettype(A,isGlobal,ierr)
    !call mat_create(W3,n1,m1+m2,k1,isGlobal,ierr)
    !call mat_zeros(W3,ierr)
    call mat_yjoin(W1,n1,m1,k1,W2,n2,m2,k2,W3,ierr)
    !print *, "W3="
    !call mat_view(W3,ierr)
    !call mat_destroy(C,ierr)
    call mat_trans(W3,C,ierr)

    call mat_destroy(W1,ierr)
    call mat_destroy(W2,ierr)
    call mat_destroy(W3,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine 

!   subroutine bk_mat_yjoin(A,m1,n1,k1,B,m2,n2,k2,C,ierr)
!     implicit none 
! #include <petsc/finclude/petscsys.h>
! #include <petsc/finclude/petscvec.h>
! #include <petsc/finclude/petscvec.h90>
! #include <petsc/finclude/petscmat.h>
!     Mat,            intent(in)  ::  A,B  
!     PetscInt,    intent(in)::  m1,n1,k1,m2,n2,k2 
!     Mat,            intent(out) ::  C
!     PetscErrorCode, intent(out) ::  ierr 
!     Mat                         ::  W1
!     Mat                         ::  I1,I2
!     PetscInt::  ista,iend
!     PetscBool    ::  isGlobal 
!     PetscScalar::  alpha
!     PetscInt:: nrow,xpos1,ypos1,xpos2,ypos2
!     integer:: i,j
!     PetscLogEvent               ::  ievent
!     call PetscLogEventRegister("mat_yjoin",0, ievent, ierr)
!     call PetscLogEventBegin(ievent,ierr)
!     ! Generate two block matrixs I1 and I2, where m1*m1 means that there are m1 rows and m1 cols in the block martrix.
!     ! I1=[I(m1*m1) 0(m1*m1)]
!     !    [0(m2*m1) 0(m2*m1)]
!     !    [0(m1*m1) I(m1*m1)]
!     !    [0(m2*m1) 0(m2*m1)]
!     !
!     ! I2=[0(m1*m2) 0(m1*m2)]
!     !    [I(m2*m2) 0(m2*m2)]
!     !    [0(m1*m2) 0(m1*m2)]
!     !    [0(m2*m2) I(m2*m2)]
!     call mat_gettype(A,isGlobal,ierr)
!     nrow=m1+m2
!     call mat_create(I1,nrow,m1,k1,isGlobal,ierr)
!     call mat_create(I2,nrow,m2,k1,isGlobal,ierr)
!     call mat_zeros(I1,ierr)
!     call mat_zeros(I2,ierr)
!     call MatGetOwnershipRange(I1,ista,iend,ierr)
!     alpha=real(1.0,kind=8)
!     do i=ista,iend-1
!        xpos1=mod(i,nrow)
!        if(xpos1<m1)  then
!           ypos1=(i/nrow)*m1+xpos1
!           call MatSetValues(I1,1,i,1,ypos1,alpha,INSERT_VALUES,ierr)
!        endif

!        xpos2=mod(i,nrow)
!        if(xpos2>=m1) then
!           ypos2=(i/nrow)*m2+mod(xpos2-m1,m2)
!           call MatSetValues(I2,1,i,1,ypos2,alpha,INSERT_VALUES,ierr)
!        endif
!     enddo
!     call mat_assemble(I1,ierr)
!     call mat_assemble(I2,ierr)
!     !print *, "I1="
!     !call mat_view(I1,ierr)
!     !print *, "I2="
!     !call mat_view(I2,ierr)

!     call mat_assemble(A,ierr)
!     call mat_assemble(B,ierr)
!     call MatMatMult(I1,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,W1,ierr)    
!     !print *, "W1="
!     !call mat_view(W1,ierr)
!     !print *, "B1="
!     !call mat_view(B,ierr)


!     call MatMatMult(I2,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,C,ierr)
!     !print *, "C1="
!     !call mat_view(C,ierr)
!     call MatAXPY(C,alpha,W1,DIFFERENT_NONZERO_PATTERN,ierr)

!     call mat_destroy(I1,ierr)
!     call mat_destroy(I2,ierr)
!     call mat_destroy(W1,ierr)
!     call PetscLogEventEnd(ievent,ierr)
!   end subroutine bk_mat_yjoin


  ! -----------------------------------------------------------------------
  !> C=[A 0] 
  !>   [0 B] 
  ! -----------------------------------------------------------------------
  subroutine mat_zjoin(A,m1,n1,k1,B,m2,n2,k2,C,ierr)
    implicit none 
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,            intent(in)  ::  A,B  
    PetscInt,  intent(in)  ::  m1,n1,k1,m2,n2,k2 
    Mat,            intent(out) ::  C
    PetscErrorCode, intent(out) ::  ierr 
    Mat                         ::  W1,W2,W3,W4
    PetscBool    ::  isGlobal 
    PetscLogEvent               ::  ievent
    call PetscLogEventRegister("mat_zjoin",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_gettype(A,isGlobal,ierr)
    call mat_create(W1,m1*k2,n1*k1,1,isGlobal,ierr)
    call mat_create(W2,m1*k1,n1*k2,1,isGlobal,ierr)
    call mat_zeros(W1,ierr)
    call mat_zeros(W2,ierr)
    !print *, "W1="
    !call mat_view(W1,ierr)
    !print *, "W2="
    !call mat_view(W2,ierr)
    !print *, "A="
    !call mat_view(A,ierr)
    !print *, "B="
    !call mat_view(B,ierr)

    call mat_xjoin(A,m1*k1,n1*k1,1,W1,m1*k2,n1*k1,1,W3,ierr)
    call mat_xjoin(W2,m1*k1,n1*k2,1,B,m1*k2,n1*k2,1,W4,ierr)
    !print *, "W3="
    !call mat_view(W3,ierr)
    !print *, "W4="
    !call mat_view(W4,ierr)
    call mat_yjoin(W3,m1*(k1+k2),n1*k1,1,W4,m1*k1,n1*k2,1,C,ierr)

    call mat_destroy(W1,ierr)
    call mat_destroy(W2,ierr)
    call mat_destroy(W3,ierr)
    call mat_destroy(W4,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine 

  ! -----------------------------------------------------------------------
  ! C=A*B
  ! -----------------------------------------------------------------------
  subroutine mat_mult(A,B,C,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::  A,B 
    Mat,intent(out)::  C
    PetscErrorCode,intent(out)::  ierr
    PetscLogEvent                   ::  ievent
    call PetscLogEventRegister("mat_mult",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(A,ierr)
    call mat_assemble(B,ierr)
    call MatMatMult(A,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,C,ierr) 

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_mult


  ! -----------------------------------------------------------------------
  ! C=A.*B
  ! -----------------------------------------------------------------------
  subroutine mat_emult(A,m1,n1,k1,B,m2,n2,k2,C,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::  A,B 
    PetscInt,intent(in)::m1,n1,k1,m2,n2,k2
    Mat,intent(out)::C
    PetscErrorCode,intent(out)::ierr
    PetscInt::col1,col2,m,n
    PetscInt,allocatable::idxn1(:),idxn2(:),idxn3(:),idxtmp(:)
    PetscScalar,allocatable::row1(:),row2(:),row3(:),rowtmp(:)
    PetscInt::  ista,iend
    PetscInt                    ::  pos1,pos2,counter
    PetscBool::  isGlobal
    integer::i
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_emult",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(A,ierr)
    call mat_assemble(B,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)
    call mat_gettype(A,isGlobal,ierr)
    call mat_create(C,m1,n1,k1,isGlobal,ierr)
    do i=ista,iend-1
       call MatGetRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
       m=col1
       call MatRestoreRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)

       call MatGetRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
       n=col2
       call MatRestoreRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)

       allocate(idxn1(m),row1(m))
       allocate(idxtmp(m),rowtmp(m))
       allocate(idxn2(n),row2(n))
       allocate(idxn3(min(m,n)),row3(min(m,n)))

       call MatGetRow(A,i,col1,idxn1,row1,ierr)
       m=col1
       idxtmp=idxn1
       rowtmp=row1
       call MatRestoreRow(A,i,col1,idxn1,PETSC_NULL_SCALAR,ierr)

       call MatGetRow(B,i,col2,idxn2,row2,ierr)
       counter=0
       pos1=1
       pos2=1
       do while(pos1<=m .and. pos2<=n)
          if(idxtmp(pos1)<idxn2(pos2)) then
             pos1=pos1+1
          elseif(idxtmp(pos1)==idxn2(pos2))then
             counter=counter+1
             idxn3(counter)=idxn2(pos2)
             row3(counter)=rowtmp(pos1)*row2(pos2)
             pos1=pos1+1
             pos2=pos2+1
          else
             pos2=pos2+1
          endif
       enddo
       call MatRestoreRow(B,i,col2,idxn2,row2,ierr)

       !print *,">i=",i,"idxn3=",idxn3,"row3=",row3
       call MatSetValues(C,1,i,counter,idxn3(1:counter),row3(1:counter),INSERT_VALUES,ierr)
       deallocate(idxn1,idxn2,idxn3,idxtmp,row1,row2,row3,rowtmp)
    enddo

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_emult


  ! -----------------------------------------------------------------------
  ! C=A./B
  ! -----------------------------------------------------------------------
  subroutine mat_ediv(A,m1,n1,k1,B,m2,n2,k2,C,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::  A,B 
    Mat,intent(out)::C
    PetscErrorCode,intent(out)::ierr
    PetscInt::col1,col2,m,n
    PetscInt::  m1,n1,k1,m2,n2,k2
    PetscInt,allocatable::idxn1(:),idxn2(:),idxn3(:),idxtmp(:)
    PetscScalar,allocatable::row1(:),row2(:),row3(:),rowtmp(:)
    PetscInt::  ista,iend
    PetscInt                    ::  pos1,pos2,counter
    PetscBool:: isGlobal
    integer::i
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_ediv",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(A,ierr)
    call mat_assemble(B,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)
    call mat_gettype(A,isGlobal,ierr)
    call mat_create(C,m1,n1,k1,isGlobal,ierr)

    do i=ista,iend-1
       call MatGetRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
       m=col1
       call MatRestoreRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)

       call MatGetRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
       n=col2
       call MatRestoreRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)

       allocate(idxn1(m),row1(m))
       allocate(idxtmp(m),rowtmp(m))
       allocate(idxn2(n),row2(n))
       allocate(idxn3(min(m,n)),row3(min(m,n)))

       call MatGetRow(A,i,col1,idxn1,row1,ierr)
       m=col1
       idxtmp=idxn1
       rowtmp=row1
       call MatRestoreRow(A,i,col1,idxn1,PETSC_NULL_SCALAR,ierr)

       call MatGetRow(B,i,col2,idxn2,row2,ierr)
       counter=0
       pos1=1
       pos2=1
       do while(pos1<=m .and. pos2<=n)
          if(idxtmp(pos1)<idxn2(pos2)) then
             pos1=pos1+1
          elseif((idxtmp(pos1)==idxn2(pos2)) .and. (row2(pos2)/=0))then
             counter=counter+1
             idxn3(counter)=idxn2(pos2)
             row3(counter)=rowtmp(pos1)/row2(pos2)
             pos1=pos1+1
             pos2=pos2+1
          else
             pos2=pos2+1
          endif
       enddo
       call MatRestoreRow(B,i,col2,idxn2,row2,ierr)

       !print *,">i=",i,"idxn3=",idxn3,"row3=",row3
       call MatSetValues(C,1,i,counter,idxn3(1:counter),row3(1:counter),INSERT_VALUES,ierr)
       deallocate(idxn1,idxn2,idxn3,idxtmp,row1,row2,row3,rowtmp)
    enddo

    call MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_ediv

  !---------------------------------------------------------
  ! replicate matrix A in x-direction for rn times, i.e.
  !  B = [A, A, ...]
  !---------------------------------------------------------
  subroutine mat_repx(A,m1,n1,k1,rn,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>

    Mat,intent(in)::  A 
    PetscInt,intent(in)::m1, n1, k1, rn
    Mat,intent(out)::B
    PetscErrorCode,intent(out)::ierr
    PetscInt::nrow,ncol
    PetscBool           ::isGlobal
    PetscInt                        :: ista, iend
    PetscInt, allocatable :: idxn1(:), idxn2(:), unbiased_idx(:)
    PetscScalar, allocatable :: row1(:), row2(:)
    PetscInt :: tmp_col, col
    PetscLogEvent :: ievent
    integer :: tmp_idx
    integer :: i, j, k

    call PetscLogEventRegister("mat_repx", 0, ievent, ierr)
    call PetscLogEventBegin(ievent, ierr)

    call MatGetSize(A,nrow,ncol,ierr)
    call mat_assemble(A, ierr)
    call MatGetOwnershipRange(A, ista, iend, ierr)
    call mat_gettype(A,isGlobal,ierr)
    call mat_create(B,  m1, n1*rn, k1, isGlobal, ierr)

    do i=ista,iend - 1

       call MatGetRow(A, i, tmp_col, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)
       col = tmp_col
       call MatRestoreRow(A, i, tmp_col, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)

       allocate(idxn1(col), row1(col), unbiased_idx(col))
       allocate(idxn2(rn*col), row2(rn*col))

       call MatGetRow(A, i, tmp_col, idxn1, row1, ierr)

       k = i / m1
       unbiased_idx = idxn1 - k * n1
       
       do j = 1,rn
          tmp_idx = (j-1)*col + 1
          idxn2(tmp_idx : tmp_idx+col-1) = k*n1*rn + (j-1)*n1 + unbiased_idx
          row2(tmp_idx  : tmp_idx+col-1) = row1
       enddo
       
       call MatRestoreRow(A, i, tmp_col, idxn1, row1, ierr)
       call MatSetValues(B, 1, i, rn*col, idxn2, row2, INSERT_VALUES, ierr)
       
       deallocate(idxn1, row1, idxn2, row2, unbiased_idx)
    enddo

    call PetscLogEventEnd(ievent, ierr)
  end subroutine mat_repx


  ! -----------------------------------------------------------------------
  ! The mat_rep function is used to replicate matrix with m times in row and
  ! n times in column. It's name refers to the repmat function in MATLAB. 
  ! Suppose P is an extended identity mM*M matrix and T is another 
  ! extended N*nN identity matrix, we can compute W=P*A firstly and then compute
  ! B=W*T. These two stpes are faster than computing B=P*A*T directly.
  ! If the size of A is M*N=3*2, suppose m=3 and n=2, we have
  ! P= [1 0 0]
  ! [0 1 0]
  ! [0 0 1]
  !    [1 0 0]
  !    [0 1 0]
  !    [0 0 1]
  !    [1 0 0]
  !    [0 1 0]
  !    [0 0 1]
  !
  ! T= [1 0 1 0]
  ! [0 1 0 1]
  ! -----------------------------------------------------------------------
  subroutine mat_rep(A,m1,n1,k1,rm,rn,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>

    Mat,intent(in)::  A 
    PetscInt,intent(in)::  m1,n1,k1,rm,rn
    Mat,intent(out)::B
    Mat :: W1,W2,W3
    PetscErrorCode,intent(out)::ierr
    PetscInt::nrow,ncol
    PetscBool           ::isGlobal
    Mat::  P,T,W
    PetscLogEvent            ::  ievent
    
    call PetscLogEventRegister("mat_rep",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_repx(A, m1, n1, k1, rn, W1, ierr)
    call mat_trans(W1, W2, ierr)
    call mat_repx(W2, rn*n1, m1, k1,rm,W3,ierr)
    call mat_trans(W3, B, ierr)

    call mat_destroy(W1, ierr)
    call mat_destroy(W2, ierr)
    call mat_destroy(W3, ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_rep

  !----------------------------------------------
  ! This function make sum along x-direction, e.g.
  !
  ! Input = [1, 2, 3]
  !         [4, 5, 6]
  ! Output = [6]
  !          [15]
  !---------------------------------------------
  subroutine mat_sumx(A, nx, ny, nz, B, ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::A
    PetscInt, intent(in) :: nx,ny,nz
    Mat,intent(out)::B
    PetscErrorCode,intent(out)::ierr
    Mat            ::W
    PetscBool      ::  isGlobal 
    PetscInt::nrow,ncol
    PetscInt::ista,iend
    PetscLogEvent            ::  ievent
    PetscInt :: col, tmp_col
    PetscScalar, allocatable :: row(:)
    PetscInt, allocatable :: idxn(:)
    PetscScalar :: sum
    integer :: i, j, k

    call PetscLogEventRegister("mat_sumx", 0, ievent, ierr)
    call PetscLogEventBegin(ievent, ierr)
    call MatGetSize(A, nrow, ncol)
    call mat_assemble(A, ierr)
    call MatGetOwnershipRange(A, ista, iend, ierr)
    call mat_gettype(A, isGlobal, ierr)
    call mat_create(B, nx, 1, nz, isGlobal, ierr)
    
    do i=ista,iend-1
       call MatGetRow(A, i, tmp_col, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)
       col = tmp_col
       call MatRestoreRow(A, i, tmp_col, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)
       allocate(idxn(col), row(col))
       call MatGetRow(A, i, tmp_col, idxn, row, ierr)

       sum = 0
       k = i / nx
       do j=1,col
          sum = sum + row(j)
       enddo
       call MatRestoreRow(A, i, tmp_col, idxn, row, ierr)
       call MatSetValues(B, 1, i, 1, k, sum, INSERT_VALUES, ierr)
       deallocate(idxn, row)
    enddo
    
    call PetscLogEventEnd(ievent, ierr)
  end subroutine mat_sumx
  
  ! -----------------------------------------------------------------------
  ! Sum of elements along with the row or column.
  ! Suppose A=[1,2,3]
  !           [4,5,6],
  ! then mat_sum(A,1,B) will make B=[5,7,9],
  !      mat_sum(A,2,B) will make B=[6 ]
  !                                 [15]
  ! -----------------------------------------------------------------------
  subroutine mat_sum(A,nx,ny,nz,ndim,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::A
    PetscInt,       intent(in)  ::  ndim
    PetscInt, intent(in) :: nx, ny, nz
    Mat,intent(out)::B
    PetscErrorCode,intent(out)::ierr
    Mat            ::W1, W2
    PetscBool      ::  isGlobal 
    PetscInt::nrow,ncol
    PetscLogEvent            ::  ievent

    call PetscLogEventRegister("mat_sum",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    if ( ndim == 1 ) then
       call mat_trans(A, W1, ierr)
       call mat_sumx(W1, ny, nx, nz, W2, ierr)
       call mat_trans(W2, B, ierr)
       
       call mat_destroy(W1,ierr)
       call mat_destroy(W2,ierr)
    else
       call mat_sumx(A, nx, ny, nz, B ,ierr)
    end if
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_sum


  ! -----------------------------------------------------------------------
  ! Compute Y = a*X + Y.
  ! -----------------------------------------------------------------------
  subroutine mat_axpy(Y,a,X,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)    ::  X 
    PetscScalar,    intent(in)    ::a
    Mat,intent(inout)   ::Y
    PetscErrorCode,intent(out)    ::ierr
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_axpy",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(X,ierr)
    call mat_assemble(Y,ierr)

    call MatAXPY(Y,a,X,DIFFERENT_NONZERO_PATTERN,ierr)
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_axpy


  ! -----------------------------------------------------------------------
  ! Compute Y = a*Y + X.
  ! -----------------------------------------------------------------------
  subroutine mat_aypx(Y,a,X,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)    ::  X 
    PetscScalar,    intent(in)    ::a
    Mat,intent(inout)   ::Y
    PetscErrorCode,intent(out)    ::ierr
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_aypx",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    call mat_assemble(X,ierr)
    call mat_assemble(Y,ierr)

    call MatAYPX(Y,a,X,DIFFERENT_NONZERO_PATTERN,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_aypx


  ! -----------------------------------------------------------------------
  ! B = A^T.
  ! -----------------------------------------------------------------------
  subroutine mat_trans(A,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)    ::  A 
    Mat,intent(out)     ::B
    PetscErrorCode,intent(out)    ::ierr
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_trans",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(A,ierr)
    call MatTranspose(A,MAT_INITIAL_MATRIX,B,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_trans


  ! -----------------------------------------------------------------------
  ! B = X*Y^T
  ! -----------------------------------------------------------------------
  subroutine mat_xyt(X,Y,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)    ::  X,Y 
    Mat,intent(out)     ::B
    Mat                 ::W
    PetscErrorCode,intent(out)    ::ierr
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_xyt",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    !MatMatTransposeMult not supported for A of type mpiaij
    !call MatMatTransposeMult(X,Y,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,B,ierr) 
    call mat_trans(Y,W,ierr)
    call mat_mult(X,W,B,ierr)
    call mat_destroy(W,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_xyt


  ! -----------------------------------------------------------------------
  ! B = X^T*Y
  ! -----------------------------------------------------------------------
  subroutine mat_xty(X,Y,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)    ::  X,Y 
    Mat,intent(out)     ::B
    PetscErrorCode,intent(out)    ::ierr
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_xty",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(X,ierr)
    call mat_assemble(Y,ierr)
    call MatTransposeMatMult(X,Y,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,B,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_xty


  ! -----------------------------------------------------------------------
  ! X = a*X
  ! -----------------------------------------------------------------------
  subroutine mat_scale(X,a,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(inout)::  X 
    PetscScalar,    intent(in)     ::a
    PetscErrorCode,intent(out)    ::ierr
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_scale",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(X,ierr)
    call MatScale(X,a,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_scale


  ! -----------------------------------------------------------------------
  ! B=fun(A,opt) 
  ! opt options: exp,log,sin,cos,tan
  ! -----------------------------------------------------------------------
  subroutine mat_math(A,opt,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_type.h"
    Mat,intent(in)::  A
    Integer,        intent(in)  ::  opt
    Mat,intent(out)::B
    PetscErrorCode,intent(out)::ierr

    PetscInt::nrow,ncol
    PetscInt::col,m
    PetscInt,allocatable        ::idxn(:),idxtmp(:)
    PetscScalar,allocatable     ::row(:),rowtmp(:)
    PetscInt::  ista,iend
    integer::i,j
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_math",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(A,ierr)
    call MatGetSize(A,nrow,ncol,ierr)

    !call mat_create(B,nrow,ncol,ierr)
    call mat_copy(A,B,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)

    do i=ista,iend-1
       call MatGetRow(A,i,col,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
       m=col
       call MatRestoreRow(A,i,col,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)

       allocate(idxn(m),idxtmp(m),row(m),rowtmp(m))

       call MatGetRow(A,i,col,idxtmp,rowtmp,ierr)
       m=col
       idxn=idxtmp
       do j=1,m
          select case(opt)
          case (MAT_MATH_EXP)
             row=exp(rowtmp)
          case (MAT_MATH_SQRT)
             row=sqrt(rowtmp)
          case (MAT_MATH_LOG)
             row=log(rowtmp)
          case (MAT_MATH_SIN)
             row=sin(rowtmp)
          case (MAT_MATH_COS)
             row=cos(rowtmp)
          case (MAT_MATH_TAN)
             row=tan(rowtmp)
          case (MAT_MATH_SQU)
             row=rowtmp**2
          case (MAT_MATH_CUBE)
             row=rowtmp**3
          case default
             row=0.0    
          end select
       enddo
       do j=1,m
          !Maybe there is a potenial bug here. When rowtmp < 1E-16, exp(rowtmp) will be NaN. 
          !To avoid this, we have to set row=0.0 at his situation.
          if(isnan(row(j))) row(j)=0.0
          !print *,"------rowtmp(",j,")=",rowtmp(j),"row(",j,")=",row(j)
       enddo

       call MatRestoreRow(A,i,col,idxtmp,rowtmp,ierr)

       !print *,">i=",i,"idxn=",idxn,"row=",row
       call MatSetValues(B,1,i,m,idxn,row,INSERT_VALUES,ierr)
       deallocate(idxn,idxtmp,row,rowtmp)
    enddo

    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_math


  ! -----------------------------------------------------------------------
  ! Solve Ax=b
  ! Dimension of A = (m, n, k)
  ! Dimension of x = (n, 1, k)
  ! Dimension of b = (m, 1, k)
  ! -----------------------------------------------------------------------
  subroutine mat_solve(A,b,x,m,n,k,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
    Mat,intent(in)::A
    Mat,intent(in)::b
    Mat,intent(out)::x
    PetscErrorCode,intent(out)::ierr
    PetscInt,intent(in)::  m,n,k
    Vec                         ::  vec_b
    Vec                         ::  vec_x
    KSP                         ::  ksp
    PC                          ::  pc
    PetscReal                   ::  tol
    PetscBool::isGlobal
    !PetscInt:: its
    PetscLogEvent            ::  ievent
    PetscInt :: nrow, ncol
    call PetscLogEventRegister("mat_solve",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    !PetscInt                    ::  its

    call mat_assemble(A,ierr)
    call mat_assemble(b,ierr)
    call mat_gettype(A,isGLobal,ierr)

    call mat_mat2vec(b, m, isGlobal,vec_b,ierr)
    
    ! if(isGlobal) then
    !    call VecCreate(PETSC_COMM_WORLD, vec_x, ierr)
    ! else
    !    call VecCreate(PETSC_COMM_SELF, vec_x, ierr)
    ! endif

    ! call VecSetSizes(vec_x, PETSC_DECIDE, n * k, ierr)
    ! call VecSetFromOptions(vec_x, ierr)

    ! call VecGetSize(vec_b, nrow, ierr)
    ! print*, "size(vec_b)=", nrow
    ! call VecGetSize(vec_x, nrow, ierr)
    ! print*, "size(vec_x)=", nrow
    ! call MatGetSize(A, nrow, ncol, ierr)
    ! print*, "size(A) = (",nrow,",",ncol,")"
    
    call VecDuplicate(vec_b,vec_x,ierr)
    
    if(isGlobal) then
       call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
       call KSPSetOperators(ksp,A,A,ierr)
       call KSPGetPC(ksp,pc,ierr)
       !call PCSetType(pc,PCBJACOBI,ierr)
    else
       call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
       call KSPSetOperators(ksp,A,A,ierr)
       call KSPGetPC(ksp,pc,ierr)
       !call PCSetType(pc,PCJACOBI,ierr)
    endif
    
    tol = 1.0e-15
    call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetFromOptions(ksp,ierr)
    call KSPSolve(ksp,vec_b,vec_x,ierr)
    
    ! call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    ! !call KSPView(ksp,PETSC_VIEWER_STDOUT_SELF,ierr)
    ! !call KSPGetIterationNumber(ksp,its,ierr)
    ! !print *, ">Iterations number=",its 
    call mat_vec2mat(vec_x, m, isGlobal, x, ierr)

    ! call mat_view(A, ierr)
    ! call sleep(1)
    ! call vec_view(vec_x, ierr)
    ! call sleep(1)
    ! call vec_view(vec_b, ierr)
    
    call KSPDestroy(ksp,ierr)
    call VecDestroy(vec_b,ierr)
    call VecDestroy(vec_x,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_solve


  ! -----------------------------------------------------------------------
  ! convert one m*1*k matrix, dimension must be (m, 1, k), into a vector m*k 
  ! -----------------------------------------------------------------------
  subroutine mat_mat2vec(A,m,isGlobal,v,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,    intent(in)      :: A 
    PetscBool,intent(in):: isGlobal
    PetscInt, intent(in) :: m
    Vec,    intent(out)     :: v
    PetscErrorCode      ::ierr

    PetscInt            ::  nrow,ncol
    PetscInt::  ista,iend
    PetscInt::  ni
    PetscInt,allocatable::  idx(:) , idy(:)
    PetscScalar,allocatable::  y(:) 
    integer :: i,pre_ik, cnt, ik

    call MatGetSize(A,nrow,ncol,ierr)

    if(isGlobal) then
       call VecCreate(PETSC_COMM_WORLD,v,ierr)
    else 
       call VecCreate(PETSC_COMM_SELF,v,ierr)
    endif
    
    call VecSetSizes(v,PETSC_DECIDE,nrow,ierr)
    call VecSetFromOptions(v,ierr)
    !call VecSetup(v,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)
    ni=iend-ista

    !print *, "ista=",ista,"iend=",iend,"ni=",ni,"==idx=",idx
    allocate(idx(ni),idy(1), y(ni))
    pre_ik = 0
    cnt = 0

    if(ista/m .ne. (iend-1)/m) then
       do i=ista,iend-1
          ik = i / m

          if((pre_ik .ne. ik) .or. (i == iend-1)) then
             print*,"pre_ik=",pre_ik, "ik=", ik
             call MatGetValues(A, cnt, idx, 1, ik, y, ierr)
             call VecSetValues(v, cnt, idx, y, INSERT_VALUES,ierr)  
             cnt = 0
             idx = 0
          endif
          
          cnt = cnt + 1          
          idx(cnt) = i
          pre_ik = ik
       enddo
    else
       ik = ista / m
       do i=ista,iend-1
          idx(i-ista+1) = i
       enddo
       idy(1) = ik
       call MatGetValues(A, ni, idx, 1, idy, y, ierr)
       call VecSetValues(v, ni, idx, y, INSERT_VALUES,ierr)  
    endif
    
    call VecAssemblyBegin(v,ierr)
    call VecAssemblyEnd(v,ierr)

    deallocate(idx, idy, y)

  end subroutine mat_mat2vec


  ! -----------------------------------------------------------------------
  ! convert one m*1 vector into a matrix 
  ! -----------------------------------------------------------------------
  subroutine mat_vec2mat(v,m, isGlobal,A,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Vec,    intent(in)      ::  v
    PetscBool, intent(in)::isGlobal
    PetscInt, intent(in) :: m
    Mat,    intent(out)     ::  A 
    PetscErrorCode              ::  ierr

    PetscInt            ::  nrow
    PetscInt::  ista,iend
    PetscInt::  ni
    PetscScalar::  y 
    integer :: i,ik,k

    call VecGetSize(v,nrow,ierr)
    k = nrow / m
    call mat_create(A,m,1,k,isGlobal,ierr)
    call mat_zeros(A,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)
    ni=iend-ista
    do i=ista,iend-1
       ik = i / m
       call VecGetValues(v,1,i,y,ierr)
       !if(y .ne. 0) then
          call MatSetValues(A,1,i,1,ik,y,INSERT_VALUES,ierr)
       !endif
    enddo
    !call mat_assemble(A,ierr)
  end subroutine mat_vec2mat

  subroutine nc_check(err, pre_msg)
    use mpi
    use pnetcdf
    implicit none
    integer err
    character(len=*) pre_msg
    ! It is a good idea to check returned value for possible error
    if (err .NE. NF90_NOERR) then
       write(6,*) trim(pre_msg), trim(nf90mpi_strerror(err))
       call MPI_Abort(MPI_COMM_WORLD, -1, err)
    end if
  end subroutine 

  !----------------------------------------------------------------
  !> Load a 3d matrix from a nc file
  !----------------------------------------------------------------
  subroutine mat_load(filename, varname, A, gnx, gny, gnz, &
       isGlobal, rank, size, ierr)
    use pnetcdf    
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname
    Mat, intent(out) :: A
    integer, intent(out) :: gnx, gny, gnz
    logical, intent(in) :: isGlobal    
    integer, intent(in) :: rank, size    
    PetscErrorCode,    intent(out)::ierr
    
    PetscScalar,allocatable  :: row(:) 
    PetscInt, allocatable    :: idxn(:)
    PetscInt:: ista,iend
    
    integer :: i,j,fid,omode,ncid
    integer :: dimid(2)

    real(kind=8),allocatable :: buf(:,:)
    character(len=1000) :: string
    PetscLogEvent       ::  ievent
    integer :: m, n, k, col
    integer :: varid,err, attrid
    integer(kind=8) :: nx,ny,nz, global_nx, global_ny,global_nz
    integer(kind=8) :: starts(2), counts(2)
    integer(kind=8) :: malloc_size=1000000, sum_size=1000000
    integer :: DD(3)
    integer :: COMM_TYPE
    
    call PetscLogEventRegister("mat_load",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    omode = NF90_NOWRITE + NF90_64BIT_OFFSET

    if(isGlobal) then
       COMM_TYPE = MPI_COMM_WORLD
    else
       COMM_TYPE = MPI_COMM_SELF
    endif
    
    ierr = nf90mpi_open(COMM_TYPE, trim(filename), omode, &
         MPI_INFO_NULL, ncid)
    
    call nc_check(ierr, 'nf90mpi_open:')

    ierr = nf90mpi_inq_varid(ncid, varname, varid)
    call nc_check(ierr, 'In nf90mpi_inq_varid varname: ')

    ierr = nf90mpi_get_att(ncid, varid, "DIM", DD)
    call nc_check(ierr, "nf90mpi_get_att: ")
    
    gnx = DD(1)
    gny = DD(2)
    gnz = DD(3)
    
    call mat_create(A, gnx, gny, gnz, isGlobal, ierr)

    call MatGetOwnershipRange(A, ista, iend, ierr)
    nx = iend - ista
    ny = gny
    
    allocate(buf(ny, nx), idxn(ny))

    starts(1) = 1
    counts(1) = ny
    starts(2) = ista+1
    counts(2) = nx
    
    ierr = nf90mpi_get_var_all(ncid, varid, buf, starts, counts)
    call nc_check(ierr, 'In nf90mpi_get_var_all: ')

    ierr = nf90mpi_close(ncid)
    call nc_check(ierr, 'In nf90mpi_close: ')
    
    do j=1,ny
       idxn(j)=j-1
    enddo
    
    do i=ista,iend-1
       k = i / gnx
       call MatSetValues(A, 1, i, ny, idxn+k*gny, buf(:,i-ista+1), &
            INSERT_VALUES, ierr)
    enddo
    
    deallocate(buf, idxn)
    
    call PetscLogEventEnd(ievent,ierr) 
  end subroutine


  ! -----------------------------------------------------------------------
  !> Save a standard row-cloumn file into a matrix 
  ! -----------------------------------------------------------------------
  subroutine mat_save(filename, varname, A, gnx, gny, gnz, isGlobal, &
       rank, size, ierr)
    use pnetcdf    
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname
    Mat,    intent(out)::A
    integer, intent(in) :: gnx, gny, gnz
    logical, intent(in) :: isGlobal
    PetscErrorCode,    intent(out)::ierr
    PetscScalar,allocatable         :: x(:,:)
    PetscScalar,allocatable         :: row(:) 
    PetscInt, allocatable    :: idxn(:)
    PetscInt    :: my_nx, my_ny, my_nz
    PetscInt    :: ista,iend
    integer :: i,j,fid,cmode,ncid
    integer :: dimid(2)
    integer, intent(in) :: rank, size
    real(kind=8),allocatable :: buf(:,:)
    character(len=1000) :: string
    PetscLogEvent       ::  ievent
    integer :: m, n, k, col
    integer :: varid,err
    integer(kind=8) :: nx, ny,nz, global_nx, global_ny,global_nz
    integer(kind=8) :: starts(2), counts(2)
    integer(kind=8) :: malloc_size=1000000, sum_size=1000000
    integer :: COMM_TYPE
    
    call PetscLogEventRegister("mat_save",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 
    
    call mat_assemble(A, ierr)
    call MatGetOwnershipRange(A, ista, iend,ierr)
    nx = iend - ista
    ny = gny
    nz = 1

    allocate(buf(ny,  nx),row(ny), idxn(ny))

    buf = 0
    do i = ista,iend-1
       k = i / gnx
       call MatGetRow(A, i, col, idxn, row, ierr)
       do j=1,col
          buf(mod(idxn(j), ny)+1, i-ista+1) = real(row(j), kind=8)
       enddo
       call MatRestoreRow(A, i, col, idxn, row, ierr)
    enddo

    if(isGlobal) then
       COMM_TYPE=MPI_COMM_WORLD
    else
       COMM_TYPE=MPI_COMM_SELF
    endif
    
    cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
    ierr = nf90mpi_create(COMM_TYPE, filename, cmode, &
         MPI_INFO_NULL, ncid)
    call nc_check(ierr, "nf90mpi_create: ")

    global_nx = gnx*gnz
    global_ny = gny
    ierr = nf90mpi_def_dim(ncid, "y",  global_ny, dimid(1))    
    ierr = nf90mpi_def_dim(ncid, "xz", global_nx, dimid(2))
    
    call nc_check(ierr, "nf90mpi_def_dim: ")
    
    ierr = nf90mpi_def_var(ncid, varname, NF90_FLOAT, dimid, varid)
    call nc_check(ierr, "nf90mpi_def_var:")

    ierr = nf90mpi_put_att(ncid, varid, "DIM", (/gnx,gny,gnz/))
    ierr = nf90mpi_enddef(ncid)
    call nc_check(ierr, "nf90mpi_enddef:")

    starts(1) = 1
    counts(1) = ny
    starts(2) = ista+1
    counts(2) = nx
    
    ierr = nf90mpi_put_var_all(ncid, varid, buf, starts, counts)
    call nc_check(ierr, "nf90mpi_put_var_all: ")
    
    ierr = nf90mpi_close(ncid)
    call nc_check(ierr, "nf90mpi_close: ")
    
    ! check if there is any PnetCDF internal malloc residue
    ierr = nf90mpi_inq_malloc_size(malloc_size)
    call nc_check(ierr, "nf90mpi_inq_malloc_size: ")
    
    deallocate(buf,row, idxn)
    call PetscLogEventEnd(ievent,ierr) 
  end subroutine

  
  
  subroutine getfilerowcol(fid, nx, ny, nz,ierr)
    implicit none
    integer, intent(in)     :: fid 
    integer, intent(out)    :: nx, ny, nz
    integer, intent(out)    :: ierr
    character(len=1)        :: cdummy
    character(len=1000)     :: string 
    integer                 :: value
    integer                 :: i,j
    nx=0
    ny=0
    nz=0
    read(fid, *) nx,ny,nz
    if(nx == 0) nx = 1
    if(ny == 0) ny = 1
    if(nz == 0) nz = 1
    ! print*, "nx=",nx,"ny=",ny,"nz=",nz
    ! read(string,*, iostat=ierr ) (nx, j=1,i)
    ! read(string,*, iostat=ierr ) (ny, j=1,i)
    ! if ( ierr == 0 ) ny = 1
    ! read(string,*, iostat=ierr ) (nz, j=1,i)
    ! if ( ierr == 0 ) nz = 1
    rewind(fid)
  end subroutine getfilerowcol

  ! -----------------------------------------------------------------------
  ! A(m,n)=value. Note that the starting point of m and n is 1. 
  ! -----------------------------------------------------------------------
  subroutine mat_setvalue(A,nx,ny,nz,m,n,k,value,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(out)::A
    PetscInt,    intent(in)::nx,ny,nz,m,n,k
    PetscScalar,    intent(in)::value
    PetscErrorCode,intent(out)::ierr
    PetscLogEvent             ::ievent
    integer :: m1,n1
    
    call PetscLogEventRegister("mat_setvalue",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    m1 = m + k * nx
    n1 = n + k * ny
    call MatSetValue(A,m1,n1,value,INSERT_VALUES,ierr)

    !call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    !call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr) 
  end subroutine mat_setvalue


  ! -----------------------------------------------------------------------
  ! B=A(rows,cols). Note that Rows and Cols should be m*1 matrix.
  ! -----------------------------------------------------------------------
  subroutine mat_getsub(A,nxa,nya,nza,idx,idy,idz,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::A
    integer,intent(in)::nxa, nya,nza,idx(:),idy(:),idz(:)
    PetscErrorCode,intent(out)::ierr
    Mat,intent(out):: B
    Vec :: V1, V2
    Mat:: W
    IS:: ISRows,ISCols
    PetscInt::  ista1,iend1, ista2, iend2
    integer:: i,k,x,ix,c1,c2
    PetscBool:: isGlobal
    PetscInt, allocatable :: row(:), col(:), vrow(:), vcol(:)
    integer ::  nxb, nyb, nzb
    integer :: nvrow
    PetscLogEvent     ::  ievent
    PetscInt :: nrow,ncol
    
    call PetscLogEventRegister("mat_getsub",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    call mat_assemble(A,ierr)
    call mat_gettype(A,isGLobal,ierr)

    nxb = size(idx)
    nyb = size(idy)
    nzb = size(idz)

    !call mat_create(B, nxb, nyb, nzb, isGlobal, ierr)
    ! call MatCreate(PETSC_COMM_WORLD, B, ierr)
    ! call MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, nxb*nzb, nyb*nzb, ierr)
    ! call MatSetFromOptions(B, ierr)
    ! call MatSetUp(B, ierr)
    !call MatGetOwnershipRange(B, ista1, iend1, ierr)
    if(isGlobal) then
       call VecCreate(PETSC_COMM_WORLD, V1, ierr)
       call VecCreate(PETSC_COMM_WORLD, V2, ierr)
    else
       call VecCreate(PETSC_COMM_SELF, V1, ierr)
       call VecCreate(PETSC_COMM_SELF, V2, ierr)
    endif
    
    call VecSetSizes(V1, PETSC_DECIDE, nxb*nzb, ierr)
    call VecSetFromOptions(V1, ierr)
    call VecSetUp(V1, ierr)
    call VecGetOwnershipRange(V1, ista1, iend1, ierr)

    call VecSetSizes(V2, PETSC_DECIDE, nyb*nzb, ierr)
    call VecSetFromOptions(V2, ierr)
    call VecSetUp(V2, ierr)
    call VecGetOwnershipRange(V2, ista2, iend2, ierr)
    
    c1 = iend1 - ista1 
    c2 = iend2 - ista2
    
    allocate(row(nxb*nzb), col(nyb*nzb), vrow(c1), vcol(c2))

    do k=0,nzb-1
       col(k*nyb+1 : k*nyb+nyb) = idy + idz(k+1) * nya
       row(k*nxb+1 : k*nxb+nxb) = idx + idz(k+1) * nxa
    enddo

    do i=ista1,iend1-1
       vrow(i-ista1+1) = row(i+1)
    enddo
    
    do i=ista2,iend2-1
       vcol(i-ista2+1) = col(i+1)
    enddo

    if(isGlobal) then
       call ISCreateGeneral(PETSC_COMM_WORLD,size(vrow),vrow, &
            PETSC_COPY_VALUES,ISRows,ierr)
       call ISCreateGeneral(PETSC_COMM_WORLD,size(vcol),vcol, &
            PETSC_COPY_VALUES,ISCols,ierr)
    else
       call ISCreateGeneral(PETSC_COMM_SELF,size(vrow),vrow, &
            PETSC_COPY_VALUES,ISRows,ierr)
       call ISCreateGeneral(PETSC_COMM_SELF,size(vcol),vcol, &
            PETSC_COPY_VALUES,ISCols,ierr)
    endif
    
    ! call MatDestroy(W1, ierr)
    ! call MatGetSubMatrix(A,ISRows,ISCols,MAT_INITIAL_MATRIX,W1,ierr)
    ! call MatView(W1, ierr)
    
    call MatGetSubMatrix(A, ISRows, ISCols, MAT_INITIAL_MATRIX, B,ierr)
    
    !call ISView(ISRows,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !call ISView(ISCols,PETSC_VIEWER_STDOUT_WORLD,ierr)

    !call MatView(W, ierr)
    !call MatView(W,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    call VecDestroy(V1,ierr)
    call VecDestroy(V2,ierr)
    call ISDestroy(ISRows,ierr)
    call ISDestroy(ISCols,ierr)
    deallocate(row, col,vrow,vcol)
    
    call PetscLogEventEnd(ievent,ierr) 
  end subroutine mat_getsub


  subroutine bk_mat_getsub(A,Rows,Cols,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::A,Rows,Cols
    PetscErrorCode,intent(out)::ierr
    Mat,intent(out):: B
    IS:: ISRows,ISCols
    PetscLogEvent               ::  ievent
    call PetscLogEventRegister("mat_getsub",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    call mat_assemble(A,ierr)
    call mat_assemble(Rows,ierr)
    call mat_assemble(Cols,ierr)

    call mat_mat2is(Rows,ISRows,ierr)
    call mat_mat2is(Cols,ISCols,ierr)
    !call ISView(ISRows,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !call ISView(ISCols,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call MatGetSubMatrix(A,ISRows,ISCols,MAT_INITIAL_MATRIX,B,ierr)
    !call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
    !call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
    call ISDestroy(ISRows,ierr)
    call ISDestroy(ISCols,ierr)
    call PetscLogEventEnd(ievent,ierr) 
  end subroutine bk_mat_getsub


  ! -----------------------------------------------------------------------
  ! convert one m*1 matrix into a IS object. 
  ! -----------------------------------------------------------------------
  subroutine mat_mat2is(A,is,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,    intent(in)      :: A 
    IS,    intent(out)     :: is 
    PetscErrorCode      ::ierr

    PetscInt            ::  nrow,ncol
    PetscInt::  ista,iend
    PetscInt::  ni
    PetscInt,allocatable::  idx(:) 
    PetscScalar,allocatable::  y(:) 
    integer :: i
    PetscBool::isGlobal 
    call MatGetSize(A,nrow,ncol,ierr)
    call mat_gettype(A,isGlobal,ierr)

    if(ncol /=1) then
       print *, "Error in mat_mat2is: the column of A should be 1"
       stop
    endif

    call MatGetOwnershipRange(A,ista,iend,ierr)
    ni=iend-ista

    allocate(idx(ni),y(ni))
    do i=ista,iend-1
       idx(i-ista+1)=i 
    enddo

    call MatGetValues(A,ni,idx,1,0,y,ierr) 
    if(isGlobal) then
       call ISCreateGeneral(PETSC_COMM_WORLD,ni,int(y),PETSC_COPY_VALUES,is,ierr)
    else
       call ISCreateGeneral(PETSC_COMM_SELF,ni,int(y),PETSC_COPY_VALUES,is,ierr)
    endif
    deallocate(idx,y)
  end subroutine mat_mat2is


  ! -----------------------------------------------------------------------
  ! Get row size and column size from A.
  ! -----------------------------------------------------------------------
  subroutine mat_getsize(A,nrow,ncol,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,        intent(in)      :: A 
    PetscInt,intent(out)     :: nrow,ncol
    PetscErrorCode          :: ierr
    PetscLogEvent               :: ievent
    call PetscLogEventRegister("mat_getsize",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    call MatGetSize(A,nrow,ncol,ierr)

    call PetscLogEventEnd(ievent,ierr) 
  end subroutine mat_getsize


  ! -----------------------------------------------------------------------
  ! Get size and column size from A.
  ! -----------------------------------------------------------------------
  subroutine mat_getownershiprange(A,ista,iend,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,        intent(in)      :: A 
    PetscInt,intent(out)     :: ista,iend 
    PetscErrorCode          :: ierr
    PetscLogEvent               :: ievent
    call PetscLogEventRegister("mat_getrange",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    call MatGetOwnershipRange(A,ista,iend,ierr)

    call PetscLogEventEnd(ievent,ierr) 
  end subroutine mat_getownershiprange


  ! -----------------------------------------------------------------------
  ! Set local array from A.
  ! -----------------------------------------------------------------------
  subroutine mat_setvalues(A,nx, ny, nz, idxm,idxn,idxk,v,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,       intent(in)      :: A
    PetscInt, intent(in) :: nx, ny, nz
    PetscInt,intent(in):: idxm(:),idxn(:),idxk(:)
    PetscScalar,intent(in):: v(:)
    PetscInt, allocatable :: idxm1(:), idxn1(:), idxn2(:)
    PetscScalar, allocatable :: vals(:)
    PetscInt:: m,n,k,ik,im,ista,iend,rid
    PetscErrorCode          :: ierr
    PetscLogEvent               :: ievent
    call PetscLogEventRegister("mat_setvalues",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    m=size(idxm)
    n=size(idxn)
    k=size(idxk)

    if(size(v) .ne. m*n*k) then
       print*, "Error: number of elements in v error."
       stop
    endif
    allocate(idxm1(m*k), idxn1(n*k), idxn2(n))

    call MatGetOwnershipRange(A, ista, iend, ierr)

    do ik=1,k
       idxm1((ik-1)*m+1:ik*m) = idxm + idxk(ik) * nx
       idxn1((ik-1)*n+1:ik*n) = idxn + idxk(ik) * ny
    enddo
    
    do im=1,m*k
       rid  = idxm1(im)
       if(rid >= ista .and. rid < iend) then
          ik = (im-1)/m
          idxn2 = idxn1(ik*n+1:(ik+1)*n)
          !print*, "rid=",rid,"idxn2=",idxn2,"v=",v((im-1)*n+1:im*n)
          call MatSetValues(A, 1, rid, &
               n,idxn2,v((im-1)*n+1:im*n),INSERT_VALUES,ierr)
       endif
    enddo

    deallocate(idxm1, idxn1, idxn2)

    call PetscLogEventEnd(ievent,ierr) 
  end subroutine 


  ! -----------------------------------------------------------------------
  ! Get local array from A.
  ! -----------------------------------------------------------------------
  subroutine mat_getvalues(A,nx,ny,nz,idxm,idxn,idxk,v,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,       intent(in)      :: A 
    PetscInt,intent(in) :: idxm(:),idxn(:),idxk(:)
    PetscInt,intent(in)  :: nx,ny,nz
    PetscInt,allocatable :: idxm1(:), idxn1(:), idxn2(:)
    PetscScalar, intent(inout):: v(:)
    PetscScalar, allocatable    :: val(:)
    PetscErrorCode              :: ierr
    PetscInt:: m,n,k,i,ik,im,ista,iend
    PetscLogEvent               :: ievent
    
    call PetscLogEventRegister("mat_getvalues",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    m=size(idxm)
    n=size(idxn)
    k=size(idxk)

    allocate(idxm1(m*k), idxn1(n*k), idxn2(n), val(n))

    !map the indicies to 2d matrix
    do ik=1,k
       idxm1((ik-1)*m+1:ik*m) = idxm + idxk(ik) * nx
       idxn1((ik-1)*n+1:ik*n) = idxn + idxk(ik) * ny
    enddo

    call mat_assemble(A,ierr)
    call MatGetOwnershipRange(A, ista, iend, ierr)

    do i=1,m*k
       im = idxm1(i)
       if(im.ge.ista .and. im<iend) then
          !ik = im / nx
          ik = (i-1)/m
          idxn2 = idxn1(ik*n+1:(ik+1)*n) 
          call MatGetValues(A,1,im, n,idxn2,val,ierr)
          v((i-1)*n+1 : i*m*n) = val
       endif
    enddo

    !call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    !call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    deallocate(idxm1, idxn1, val)
    call PetscLogEventEnd(ievent,ierr) 
  end subroutine 



  ! -----------------------------------------------------------------------
  ! Norm(A)
  ! -----------------------------------------------------------------------
  subroutine mat_norm(A,ntype,res,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,       intent(in)      :: A
    integer,intent(in):: ntype 
    PetscReal,  intent(out):: res 
    PetscErrorCode          :: ierr
    PetscLogEvent               :: ievent
    call PetscLogEventRegister("mat_norm",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    call mat_assemble(A,ierr)
    call MatNorm(A,ntype,res,ierr) 
    call PetscLogEventEnd(ievent,ierr) 
  end subroutine mat_norm


  ! -----------------------------------------------------------------------
  ! C= A<B, A>B, A<=B, or A>=B 
  ! -----------------------------------------------------------------------
  subroutine mat_compare(A,B,opt,C,nx,ny,nz,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,        intent(in)::  A,B
    PetscInt, intent(in) :: nx, ny, nz
    integer,        intent(in)  ::  opt
    Mat,            intent(out)::C
    PetscErrorCode,intent(out)::ierr
    Mat                         ::W
    PetscInt::nrow1,ncol1,nrow2,ncol2
    PetscInt::col,coltmp
    PetscInt,allocatable::idxn(:),idxntmp(:)
    PetscScalar,allocatable::row(:),rowtmp(:)
    PetscInt::  i,ista,iend
    PetscScalar::  alpha 
    PetscBool             ::isGlobal
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_compare",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(A,ierr)
    call mat_assemble(B,ierr)
    call mat_gettype(A,isGlobal,ierr)

    call mat_create(C, nx, ny, nz, isGlobal,ierr)
    call mat_zeros(C,ierr)    
    call mat_assemble(C, ierr)
    
    alpha=-1.0
    call MatDuplicate(A,MAT_COPY_VALUES,W,ierr)
    call MatAXPY(W,alpha, B, DIFFERENT_NONZERO_PATTERN,ierr)
    call MatGetOwnershipRange(W,ista,iend,ierr)
    
    allocate(idxn(ny),row(ny), idxntmp(ny), rowtmp(ny))

    do i=ista,iend-1
       call MatGetRow(W, i, col, idxn, row,ierr)
       select case(opt)
       case (MAT_COMPARE_LT)
          where(row <  0) 
             row=1
          else where
             row=0.0
          end where
       case (MAT_COMPARE_LE)
          where(row <= 0) 
             row=1
          else where
             row=0.0
          end where
       case (MAT_COMPARE_GT)
          where(row >  0) 
             row=1
          else where
             row=0.0
          end where
       case (MAT_COMPARE_GE)
          where(row >= 0) 
             row=1
          else where
             row=0.0
          end where
       case (MAT_COMPARE_EQ)
          where(row == 0) 
             row=1
          else where
             row=0.0
          end where
       case (MAT_COMPARE_NQ)
          where(row /= 0) 
             row=1
          else where
             row=0.0
          end where
       case default
          row=0.0    
       end select
       !print*,"i=", i, "col=", col, "idxn=", idxn, "row=",row
       call MatSetValues(C, 1, i, col, idxn, row, INSERT_VALUES, ierr);
       call MatRestoreRow(W, i, col, idxn, row, ierr)
    enddo

    deallocate(idxn, row, idxntmp, rowtmp)
    !call mat_assemble(C,ierr)
    call mat_destroy(W,ierr) 
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_compare


  ! -----------------------------------------------------------------------
  ! Create sparse matrix
  ! -----------------------------------------------------------------------
  subroutine mat_sparse(Ind_m,Ind_n,Ind_k,A,m,n,k,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::  Ind_m, Ind_n, Ind_k 
    Mat,intent(in)::  A
    PetscInt,intent(in)::  m,n,k
    Mat,intent(out)::  B
    PetscErrorCode,intent(out)::  ierr
    PetscInt::  nrow_m,nrow_n,nrow_k,nrow_a,nrow
    PetscInt::  ncol_m,ncol_n,ncol_k,ncol_a,ncol
    PetscScalar::  row_m,row_n,row_k,row_a
    PetscInt::  ista,iend
    integer::  i
    PetscBool::  isGlobal
    PetscLogEvent                ::  ievent
    call PetscLogEventRegister("mat_sparse",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(Ind_m,ierr) 
    call mat_assemble(Ind_n,ierr)
    call mat_assemble(Ind_k,ierr)
    
    call mat_assemble(A,ierr) 
    call MatGetSize(Ind_m,nrow_m,ncol_m,ierr)
    call MatGetSize(Ind_n,nrow_n,ncol_n,ierr)
    call MatGetSize(Ind_k,nrow_k,ncol_k,ierr)
    call MatGetSize(A,nrow_a,ncol_a,ierr)
    
    !call mat_view(Ind_m,ierr)
    !call mat_view(Ind_n,ierr)
    !call mat_view(A,ierr)
    if(nrow_m/=nrow_n .or. nrow_m/=nrow_k .or. nrow_n/=nrow_k) then
       print *, "Error in mat_sparse: matrix Ind_m, matrix Ind_n &
            and matrix A should have the same row size"
       stop
    endif

    if(ncol_m/=1.or. ncol_m/=1 .or. ncol_m/=1 .or. ncol_k/=1) then
       print *, "Error in mat_sparse: matrix Ind_m, matrix Ind_n &
            and matrix A should have only one column"
       stop
    endif
    
    call mat_gettype(A,isGlobal,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)

    do i=ista,iend-1
       call MatGetRow(Ind_m,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row_m,ierr)
       call MatRestoreRow(Ind_m,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row_m,ierr)

       call MatGetRow(Ind_n,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row_n,ierr)
       call MatRestoreRow(Ind_n,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row_n,ierr)

       call MatGetRow(Ind_k,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row_k,ierr)
       call MatRestoreRow(Ind_k,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row_k,ierr)

       call MatGetRow(A, i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row_a,ierr)
       call MatRestoreRow(A,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row_a,ierr)
       nrow = row_m + row_k * m
       ncol = row_n + row_k * n
       !print *,">row1=",row1,"row2=",row2,"row3=",row3
       call MatSetValues(B,1,nrow,1,ncol,row_a, INSERT_VALUES, ierr)
    enddo

    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_sparse




  ! -----------------------------------------------------------------------
  !> Transform Cartesian to spherical coordinates.
  !! [TH,PHI,R] = cart2sph(X,Y,Z) transforms corresponding elements of
  !!    data stored in Cartesian coordinates X,Y,Z to spherical
  !!    coordinates (azimuth TH, elevation PHI, and radius R).
  !!    TH and PHI are returned in radians.
  !! where, \n
  !!   azimuth = atan2(y,x)
  !!   elevation = atan2(z,sqrt(x.^2 + y.^2))
  !!   r = sqrt(x.^2 + y.^2 + z.^2)
  ! -----------------------------------------------------------------------
  subroutine mat_cart2sph(A,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in):: A
    Mat,intent(out):: B
    PetscErrorCode,intent(out)        :: ierr

    PetscInt::  nrow1,ncol1,nrow2,ncol2
    PetscInt,allocatable        ::  idx1(:),idx2(:)
    PetscScalar,allocatable     ::  row1(:),row2(:)
    PetscInt::  ista,iend
    integer::  i,j
    call mat_assemble(A,ierr) 
    call MatGetSize(A,nrow1,ncol1,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)

    nrow2=nrow1
    ncol2=3

    allocate(idx1(ncol1),idx2(ncol2),row1(ncol1),row2(ncol2))

    do i=ista,iend-1
       do j=1,ncol1
          idx1(j)=j-1
       enddo
       do j=1,ncol2
          idx2(j)=j-1
       enddo

       call MatGetRow(A,i,ncol1,idx1,row1,ierr)
       row2(1) = atan2(row1(2),row1(1)) 
       row2(2) = atan2(row1(3),sqrt(row1(1)**2+row1(2)**2))
       row2(3) = sqrt(row1(1)**2+row1(2)**2+row1(3)**2) 
       call MatRestoreRow(A,i,ncol1,idx1,row1,ierr)

       call MatSetValues(B,1,i,ncol2,idx2,row2,INSERT_VALUES,ierr)
    enddo

    deallocate(idx1,idx2,row1,row2)
    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
  end subroutine mat_cart2sph


  ! -----------------------------------------------------------------------
  !> Set the diagonal of A to constant value 
  ! -----------------------------------------------------------------------
  subroutine mat_setdiag(A,value,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(inout)::A
    PetscScalar,    intent(in)::value
    PetscErrorCode,intent(out)::ierr
    Vec::  x 
    PetscInt::  nrow,ncol 
    PetscBool:: isGlobal
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_setdiag",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_gettype(A,isGlobal,ierr)
    call MatGetSize(A,nrow,ncol,ierr)
    if(nrow /= ncol) then
       print *, "Error in mat_diag_set: the row number should equal to the column number"
       stop
    endif
    if(isGlobal) then
       call VecCreate(PETSC_COMM_WORLD,x,ierr)
    else
       call VecCreate(PETSC_COMM_SELF,x,ierr)
    endif
    call VecSetSizes(x,PETSC_DECIDE,nrow,ierr)
    call VecSetFromOptions(x,ierr)
    call VecSet(x,value,ierr)
    !call VecAssemblyBegin(x,ierr)
    !call VecAssemblyEnd(x,ierr)
    !call mat_assemble(A,ierr)

    call MatDiagonalSet(A,x,INSERT_VALUES,ierr)

    call VecDestroy(x,ierr)

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_setdiag


  subroutine mat_setcol(A,idxn,B,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(inout)::  A 
    Mat,intent(in)::  B 
    PetscInt,intent(in)::  idxn
    PetscErrorCode,intent(out)::  ierr
    PetscInt::  nrow1,ncol1,nrow2,ncol2
    PetscInt::  col2,m
    PetscInt::  idxn1(1),idxn2(1)
    PetscScalar::  row1(1),row2(1)
    PetscInt::  ista,iend
    integer::  i
    PetscLogEvent            ::  ievent
    call PetscLogEventRegister("mat_setcol",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_assemble(A,ierr)
    call mat_assemble(B,ierr)

    call MatGetSize(A,nrow1,ncol1,ierr)
    call MatGetSize(B,nrow2,ncol2,ierr)
    if(nrow1/=nrow2)then
       print *, "Error in mat_setcol: Matrix A and Matrix B &
            should have the same row number"
       stop
    endif
    if(ncol2/=1)then
       print *, "Error in mat_setcol: the column number of Matrix B should be 1"
       stop
    endif

    call MatGetOwnershipRange(A,ista,iend,ierr)

    m=0
    !call MatZeroRows(A,1,idxn,0.0,0,0,ierr) 
    do i=ista,iend-1
       idxn1=0
       idxn2=0
       row2=0
       row1=0 
       call MatGetRow(B,i,col2,idxn2,row2,ierr)
       !print *,">i=",i,"idxn2=",idxn1,"row2=",row2
       m=col2
       idxn1=idxn2
       row1=row2
       call MatRestoreRow(B,i,col2,idxn2,row2,ierr)

       !print *,">m=",m,"idxn1=",idxn1,"row1=",row1
       !if(m > 0) then
       call MatSetValues(A,1,i,1,idxn,row1,INSERT_VALUES,ierr)
       !endif
    enddo

    call PetscLogEventEnd(ievent,ierr)
  end subroutine mat_setcol

  subroutine mat_gettype(A,isGlobal,ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,intent(in)::A
    PetscBool,intent(out)::isGlobal
    PetscErrorCode,intent(out)::ierr
    MatType                     ::  flag 

    call MatGetType(A,flag,ierr)  
    if(flag/=MATSEQAIJ)then   
       isGlobal=.true.
    else 
       isGlobal=.false.
    endif
  end subroutine mat_gettype


end module dm_mat
