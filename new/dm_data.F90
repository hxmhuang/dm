! -----------------------------------------------------------------------
!> Distributed Matrix Computing Library
! -----------------------------------------------------------------------
module dm_data
  use dm_common
  implicit none
#include "dm_type.h"
  
  public :: dm_init, data_plus
  
  interface get_local_ptr
     module procedure get_local_ptr1d
     module procedure get_local_ptr2d
     module procedure get_local_ptr3d
  end interface get_local_ptr

  interface data_plus
     module procedure data_plus1
     module procedure data_plus2
  end interface data_plus

  interface data_minus
     module procedure data_minus1
     module procedure data_minus2
  end interface

  interface data_mult
     module procedure data_mult1
     module procedure data_mult2
  end interface data_mult

  interface data_divd
     module procedure data_divd1
     module procedure data_divd2
  end interface data_divd

  interface data_eop
     module procedure data_eop1
     module procedure data_eop2
  end interface data_eop
  
contains

  subroutine dm_init(ierr)
#include "petsc.h"
    integer, intent(out) :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr) 
  end subroutine dm_init

  subroutine dm_finalize(ierr)
    integer, intent(out) :: ierr
    call PetscFinalize(ierr)
  end subroutine 

!   subroutine dm_clone_dm(dst, src)
!     implicit none
! #include "petsc.h"
!     DM, intent(out) :: dst
!     DM, intent(in) :: src
!     integer :: ierr
!     call DMClone(src, dst, ierr)
!   end subroutine

  subroutine data_duplicate(dst, src, ierr)
    implicit none
#include "petsc.h"
    Vec, intent(out) :: dst
    Vec, intent(in)  :: src
    ! DM, intent(out)  :: dst_dm
    ! DM, intent(in)   :: src_dm
    integer, intent(out) :: ierr
    
    !call DMClone(src_dm, dst_dm, ierr)
    call VecDuplicate(src, dst, ierr)

    ! call VecView(src,PETSC_VIEWER_STDOUT_WORLD,ierr)
    ! call VecView(dst,PETSC_VIEWER_STDOUT_WORLD,ierr)
  end subroutine 
  
  subroutine data_clone(dst, src, ierr)
    implicit none
#include "petsc.h"
    Vec, intent(out) :: dst
    Vec, intent(in)  :: src
    ! DM, intent(out)  :: dst_dm
    ! DM, intent(in)   :: src_dm
    integer, intent(out) :: ierr
    
    !call DMClone(src_dm, dst_dm, ierr)
    call VecDuplicate(src, dst, ierr)
    call VecCopy(src, dst, ierr)
    ! call VecView(src,PETSC_VIEWER_STDOUT_WORLD,ierr)
    ! call VecView(dst,PETSC_VIEWER_STDOUT_WORLD,ierr)
  end subroutine 

  subroutine data_destroy(data, ierr)
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: data
    DM :: data_dm
    integer, intent(out) :: ierr

    if(data /= 0) then
       call VecGetDM(data, data_dm, ierr)       
       call DMRestoreGlobalVector(data_dm, data, ierr)
       !call DMDestroy(data_dm, ierr)
       data = 0
    endif
  end subroutine 
  
!   subroutine dm_destroy_dm(ada)
!     implicit none
! #include "petsc.h"
!     DM, intent(inout) :: ada
!     integer :: ierr
!     call DMDestroy(ada, ierr)
!   end subroutine 

  
  subroutine get_local_ptr1d(A, xx)
    implicit none
#include "petsc.h"
    Vec, intent(in) :: A
    DM :: ada
    PetscScalar, pointer, intent(out) :: xx(:)
    integer :: ierr
    PetscScalar, pointer :: raw_ptr(:)
    PetscInt :: xs, xl

    call VecGetDM(A, ada, ierr)
    call DMDAVecGetArrayF90(ada, A, raw_ptr, ierr)
    call DMDAGetCorners(ada,xs, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         xl,PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
    xx => raw_ptr(xs:xs+xl-1)
  end subroutine

  subroutine get_local_ptr2d(A, xx)
    implicit none
#include "petsc.h"
    Vec, intent(in) :: A
    DM :: ada
    PetscScalar, pointer, intent(out) :: xx(:,:)
    integer :: ierr
    PetscScalar, pointer :: raw_ptr(:,:)
    PetscInt :: xs, xl, ys, yl

    call VecGetDM(A, ada, ierr)
    call DMDAVecGetArrayF90(ada, A, raw_ptr, ierr)
    call DMDAGetCorners(ada,xs, ys, PETSC_NULL_INTEGER, &
         xl, yl, PETSC_NULL_INTEGER, ierr)
    
    xx => raw_ptr(xs:xs+xl-1, ys:ys+yl-1)

  end subroutine

  subroutine print3d(A, prefix)
    implicit none
#include "petsc.h"
    Vec :: A, vout
    character(len=*) :: prefix
    VecScatter :: ctx
    integer :: ierr
    
    call  VecScatterCreateToZero(A, ctx, vout, ierr);
    ! scatter as many times as you need
    call  VecScatterBegin(ctx, A, vout,INSERT_VALUES,SCATTER_FORWARD, ierr)
    call  VecScatterEnd(ctx, A,vout,INSERT_VALUES,SCATTER_FORWARD, ierr)
    ! destroy scatter context and local vector when no longer needed
    call  VecScatterDestroy(ctx, ierr);
    
    call  VecDestroy(vout, ierr);

  end subroutine
  
  subroutine get_local_ptr3d(A, xx)
    implicit none
#include "petsc.h"
    Vec, intent(in) :: A
    DM :: ada
    PetscScalar, pointer, intent(out) :: xx(:,:,:)
    integer :: ierr
    PetscScalar, pointer :: raw_ptr(:,:,:)
    PetscInt :: xs, xl, ys, yl, zs, zl
    
    call VecGetDM(A, ada, ierr)
    
    call DMDAVecGetArrayF90(ada, A, raw_ptr, ierr)
    
    call DMDAGetCorners(ada,xs, ys, zs, &
         xl, yl, zl, ierr)
    
    xx => raw_ptr(xs:xs+xl-1, ys:ys+yl-1, zs:zs+zl-1)
    
    ! call VecView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    ! print*, "raw_ptr:", raw_ptr(xs:xs+xl-1, ys:ys+yl-1, zs:zs+zl-1)

  end subroutine

!   subroutine get_corners(A, cor)
!     implicit none
! #include "petsc.h"
!     Vec :: A
!     type(corner), intent(out) :: cor 
!     PetscInt :: xs, xl, ys, yl, zs, zl
!     DM :: adm
!     integer :: ierr

!     call assert(size(starts) == size(ends) &
!          .and. size(starts) > 0, &
!          __FILE__, __LINE__)
    
!     call VecGetDM(A, adm, ierr)
    
!     select case(size(starts))
!     case (1)
!        call DMDAGetCorners(adm, xs, &
!             PETSC_NULL_INTEGER, &
!             PETSC_NULL_INTEGER, &
!             xl, PETSC_NULL_INTEGER, &
!             PETSC_NULL_INTEGER, ierr)
!        cor%xs = xs
!        cor%xe = xs + xl - 1;
!     case (2)
!        call DMDAGetCorners(adm, xs, &
!             ys, &
!             PETSC_NULL_INTEGER, &
!             xl, yl, &
!             PETSC_NULL_INTEGER, ierr)
!        cor%xs = xs
!        cor%ys = ys
!        cor%xe = xs + xl - 1
!        cor%ye = ys + yl - 1
!     case (3)
!        call DMDAGetCorners(adm, starts(1), &
!             starts(2), &
!             starts(3), &
!             ends(1), &
!             ends(2), &
!             ends(3), ierr)

!     end select
!   end subroutine
     
  ! -----------------------------------------------------------------------
  ! A=1. 
  ! -----------------------------------------------------------------------
  subroutine data_constants(A, alpha, t_shape, ierr)
    implicit none
#include "petsc.h"

    Vec, intent(out):: A
    DM :: ada
    
    PetscErrorCode, intent(out) ::ierr
    integer :: t_shape(:)

    PetscInt :: m, n, k

    integer :: num_dim
    PetscInt :: xs, ys, zs, xl, yl, zl
    PetscScalar,intent(in)::  alpha
    PetscInt::  ista,iend,ilocal
    PetscInt,allocatable::idxm(:),idxn(:)
    PetscScalar,allocatable::row(:)
    integer :: i,j
    PetscLogEvent            ::  ievent
    integer :: dof = 1, s = 1

    PetscScalar, pointer :: x1(:), x2(:,:), x3(:,:,:)

    call PetscLogEventRegister("data_constants",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    num_dim = size(t_shape)

    m = t_shape(1)
    if(num_dim >= 2) n = t_shape(2)
    if(num_dim >= 3) k = t_shape(3)

    if(num_dim == 1) then
       call DMDACreate1d(PETSC_COMM_WORLD, &
         &     DM_BOUNDARY_NONE,        &
         &     m,  dof, s,     &
         &     PETSC_NULL_INTEGER,ada,ierr)
       
       call DMGetGlobalVector(ada, A,ierr)
       call DMDAGetCorners(ada,xs, PETSC_NULL_INTEGER, &
            PETSC_NULL_INTEGER, xl, &
            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)       
       call DMDAVecGetArrayF90(ada, A, x1,ierr)

       !print*, "xs=", xs, "xl=", xl
       x1(xs:xs+xl-1) = alpha;
       
       call DMDAVecRestoreArrayF90(ada, A, x1, ierr)
       !call VecView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
       !call DMRestoreGlobalVector(ada, A, ierr)
       !call DMDestroy(ada,ierr)
       
    else if (num_dim == 2) then

       call DMDACreate2d(PETSC_COMM_WORLD, &
         &     DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,        &
         &     DMDA_STENCIL_STAR, m,n, PETSC_DECIDE,      &
         &     PETSC_DECIDE,dof,s,                       &
         &     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,    &
         &     ada,ierr)
       
       call DMGetGlobalVector(ada, A,ierr)
       call DMDAGetCorners(ada,xs, ys, &
            PETSC_NULL_INTEGER, xl, yl, &
            PETSC_NULL_INTEGER, ierr)       
       call DMDAVecGetArrayF90(ada, A, x2,ierr)
       
       
       x2(xs:xs+xl-1, ys:ys+yl-1) = alpha

       call DMDAVecRestoreArrayF90(ada, A, x2, ierr)
       !call VecView(A, PETSC_VIEWER_STDOUT_WORLD,ierr)
       !call DMRestoreGlobalVector(ada, A, ierr)
       !call DMDestroy(ada,ierr)

    else if(num_dim == 3) then
       
       call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,               &
            &     DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,                      &
            &     DMDA_STENCIL_STAR, m,n,k,PETSC_DECIDE,PETSC_DECIDE,     &
            &     PETSC_DECIDE,dof,s,                                     &
            &     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                  &
            &     PETSC_NULL_INTEGER,ada,ierr)

       call DMGetGlobalVector(ada, A,ierr)
       call DMDAGetCorners(ada,xs,ys,zs, xl,yl,zl,ierr)
       call DMDAVecGetArrayF90(ada, A, x3,ierr)

       x3(zs:xs+xl-1, ys:ys+yl-1, zs:zs+zl-1) = alpha;
       
       call DMDAVecRestoreArrayF90(ada, A, x3, ierr)
       !call VecView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
       !call DMRestoreGlobalVector(ada, A, ierr)
       !call DMDestroy(ada,ierr)
    endif

    call PetscLogEventEnd(ievent,ierr)
  end subroutine data_constants

  subroutine data_rand(A, t_shape, ierr)
    implicit none
#include "petsc.h"

    Vec, intent(out):: A
    DM :: ada
    PetscErrorCode, intent(out) ::ierr
    integer :: t_shape(:)
    PetscInt :: m, n, k
    integer :: num_dim
    PetscInt :: xs, ys, zs, xl, yl, zl
    PetscInt::  ista,iend,ilocal
    PetscInt,   allocatable::idxm(:),idxn(:)
    PetscScalar,allocatable::row(:)
    integer :: i,j
    PetscLogEvent            ::  ievent
    integer :: dof = 1, s = 1
    PetscRandom :: rctx
    
    call PetscLogEventRegister("data_rand",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    num_dim = size(t_shape)

    m = t_shape(1)
    if(num_dim >= 2) n = t_shape(2)
    if(num_dim >= 3) k = t_shape(3)

    call PetscRandomCreate(PETSC_COMM_WORLD, rctx, ierr)

    select case(num_dim)
    case (1)
       call DMDACreate1d(PETSC_COMM_WORLD, &
            &     DM_BOUNDARY_NONE,        &
            &     m,  dof, s,     &
            &     PETSC_NULL_INTEGER,ada,ierr)
    case (2)
       call DMDACreate2d(PETSC_COMM_WORLD, &
            &     DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,        &
            &     DMDA_STENCIL_STAR, m,n, PETSC_DECIDE,     &
            &     PETSC_DECIDE,dof,s,                       &
            &     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,    &
            &     ada,ierr)

    case (3)
       call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,               &
            &     DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,                      &
            &     DMDA_STENCIL_STAR, m,n,k,PETSC_DECIDE,PETSC_DECIDE,     &
            &     PETSC_DECIDE,dof,s,                                     &
            &     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                  &
            &     PETSC_NULL_INTEGER,ada,ierr)
    end select

   call DMCreateGlobalVector(ada, A, ierr)
   call VecSetRandom(A, rctx)
   call PetscLogEventEnd(ievent,ierr)
  end subroutine 
  
  subroutine data_plus1(A, B, ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B
    integer, intent(out) :: ierr    
    call data_eop(A, B, type_plus, ierr)
  end subroutine 

  subroutine data_plus2(A, B, C,  ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B, C
    integer, intent(out) :: ierr
    call data_eop(A, B, C, type_plus, ierr)
  end subroutine 

  subroutine data_minus1(A, B, ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B
    integer, intent(out) :: ierr    
    call data_eop(A, B, type_minus, ierr)
  end subroutine 

  subroutine data_minus2(A, B, C,  ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B, C
    integer, intent(out) :: ierr
    call data_eop(A, B, C, type_minus, ierr)
  end subroutine 

  subroutine data_mult1(A, B, ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B
    integer, intent(out) :: ierr    
    call data_eop(A, B, type_mult, ierr)
  end subroutine 

  subroutine data_mult2(A, B, C,  ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B, C
    integer, intent(out) :: ierr
    call data_eop(A, B, C, type_mult, ierr)
  end subroutine 

  subroutine data_divd1(A, B, ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B
    integer, intent(out) :: ierr    
    call data_eop(A, B, type_divd, ierr)
  end subroutine 

  subroutine data_divd2(A, B, C,  ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B, C
    integer, intent(out) :: ierr
    call data_eop(A, B, C, type_divd, ierr)
  end subroutine 

  !element-wise bindary operation
  !A = A + B, A = A - B, A = A .* B, A = A ./ B
  subroutine data_eop1(A, B, op_type, ierr) 
    implicit none
#include "petsc.h"

    Vec, intent(inout) :: A
    Vec, intent(in) :: B
    DM :: A_dm, B_dm
    PetscErrorCode, intent(out) ::ierr
    PetscInt :: m, n, k
    integer  :: num_dim
    PetscLogEvent            ::  ievent
    PetscInt :: xs, xl, ys, yl, zs, zl, xe,ye,ze
    PetscScalar, pointer :: a_1d(:),     b_1d(:)
    PetscScalar, pointer :: a_2d(:,:),   b_2d(:,:)
    PetscScalar, pointer :: a_3d(:,:,:), b_3d(:,:,:)
    integer, intent(in) :: op_type
    
    call PetscLogEventRegister("data_ele_bop",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    
    call VecGetDM(A, A_dm, ierr)
    call VecGetDM(A, B_dm, ierr)

    call DMDAGetCorners(A_dm, xs, ys, zs, xl, yl, zl,ierr)
    xe = xs + xl - 1
    ye = ys + yl - 1
    ze = zs + zl - 1
    
    call DMDAGetInfo(A_dm, num_dim, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

    select case (num_dim)
    case (1)
       call DMDAVecGetArrayF90(A_dm, A, a_1d, ierr)
       call DMDAVecGetArrayF90(B_dm, B, b_1d, ierr)

       select case(op_type)
       case (type_plus)
          a_1d(xs:xe) = a_1d(xs:xe) + b_1d(xs:xe)
       case (type_minus)
          a_1d(xs:xe) = a_1d(xs:xe) - b_1d(xs:xe)
       case (type_mult)
          a_1d(xs:xe) = a_1d(xs:xe) * b_1d(xs:xe)
       case (type_divd)
          a_1d(xs:xe) = a_1d(xs:xe) / b_1d(xs:xe)
       end select

       call DMDAVecRestoreArrayF90(A_dm, A, a_1d, ierr)
       call DMDAVecRestoreArrayF90(B_dm, B, b_1d, ierr)
    case (2)
       call DMDAVecGetArrayF90(A_dm, A, a_2d, ierr)
       call DMDAVecGetArrayF90(B_dm, B, b_2d, ierr)

       select case(op_type)
       case (type_plus)
          a_2d(xs:xe, ys:ye) = a_2d(xs:xe, ys:ye) + b_2d(xs:xe, ys:ye)
       case (type_minus)
          a_2d(xs:xe, ys:ye) = a_2d(xs:xe, ys:ye) - b_2d(xs:xe, ys:ye)
       case (type_mult)
          a_2d(xs:xe, ys:ye) = a_2d(xs:xe, ys:ye) * b_2d(xs:xe, ys:ye)
       case (type_divd)
          a_2d(xs:xe, ys:ye) = a_2d(xs:xe, ys:ye) / b_2d(xs:xe, ys:ye)
       end select

       call DMDAVecRestoreArrayF90(A_dm, A, a_2d, ierr)
       call DMDAVecRestoreArrayF90(B_dm, B, b_2d, ierr)
    case (3)
       call DMDAVecGetArrayF90(A_dm, A, a_3d, ierr)
       call DMDAVecGetArrayF90(B_dm, B, b_3d, ierr)

       select case(op_type)
       case (type_plus)
          a_3d(xs:xe, ys:ye, zs:ze) = a_3d(xs:xe, ys:ye, zs:ze) + &
               b_3d(xs:xe, ys:ye, zs:ze)
       case (type_minus)
          a_3d(xs:xe, ys:ye, zs:ze) = a_3d(xs:xe, ys:ye, zs:ze) - &
               b_3d(xs:xe, ys:ye, zs:ze)
       case (type_mult)
          a_3d(xs:xe, ys:ye, zs:ze) = a_3d(xs:xe, ys:ye, zs:ze) * &
               b_3d(xs:xe, ys:ye, zs:ze)
       case (type_divd)
          a_3d(xs:xe, ys:ye, zs:ze) = a_3d(xs:xe, ys:ye, zs:ze) / &
               b_3d(xs:xe, ys:ye, zs:ze)
       end select

       call DMDAVecRestoreArrayF90(A_dm, A, a_3d, ierr)
       call DMDAVecRestoreArrayF90(B_dm, B, b_3d, ierr)

    end select
    
    call PetscLogEventEnd(ievent, ierr)
  end subroutine 

  !A = B + C
  subroutine data_eop2(A, B, C, op_type, ierr)
    implicit none
#include "petsc.h"

    Vec, intent(inout) :: A
    Vec, intent(in) :: B
    Vec, intent(in) :: C
    DM :: A_dm, B_dm, C_dm
    PetscErrorCode, intent(out) ::ierr
    PetscInt :: m, n, k
    integer  :: num_dim
    PetscLogEvent            ::  ievent
    PetscInt :: xs, xl, ys, yl, zs, zl, xe,ye,ze
    integer, intent(in) :: op_type
    PetscScalar, pointer :: a_1d(:),     b_1d(:),     c_1d(:)    
    PetscScalar, pointer :: a_2d(:,:),   b_2d(:,:),   c_2d(:,:)    
    PetscScalar, pointer :: a_3d(:,:,:), b_3d(:,:,:), c_3d(:,:,:)
    
    call PetscLogEventRegister("data_plus",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    
    call VecGetDM(A, A_dm, ierr)
    call VecGetDM(A, B_dm, ierr)
    call VecGetDM(A, C_dm, ierr)

    call DMDAGetCorners(A_dm, xs, ys, zs, xl, yl, zl,ierr)
    xe = xs + xl - 1
    ye = ys + yl - 1
    ze = zs + zl - 1
    
    call DMDAGetInfo(A_dm, num_dim, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

    select case (num_dim)
    case (1)
       call DMDAVecGetArrayF90(A_dm, A, a_1d, ierr)
       call DMDAVecGetArrayF90(B_dm, B, b_1d, ierr)
       call DMDAVecGetArrayF90(C_dm, C, c_1d, ierr)

       select case(op_type)
       case (type_plus)
          a_1d(xs:xe) = b_1d(xs:xe) + c_1d(xs:xe)
       case (type_minus)
          a_1d(xs:xe) = b_1d(xs:xe) - c_1d(xs:xe)
       case (type_mult)
          a_1d(xs:xe) = b_1d(xs:xe) * c_1d(xs:xe)
       case (type_divd)
          a_1d(xs:xe) = b_1d(xs:xe) / c_1d(xs:xe)
       end select

       call DMDAVecRestoreArrayF90(A_dm, A, a_1d, ierr)
       call DMDAVecRestoreArrayF90(B_dm, B, b_1d, ierr)
       call DMDAVecRestoreArrayF90(C_dm, C, c_1d, ierr)
    case (2)
       call DMDAVecGetArrayF90(A_dm, A, a_2d, ierr)
       call DMDAVecGetArrayF90(B_dm, B, b_2d, ierr)
       call DMDAVecGetArrayF90(C_dm, C, c_2d, ierr)

       select case(op_type)
       case (type_plus)
          a_2d(xs:xe, ys:ye) = b_2d(xs:xe, ys:ye) + c_2d(xs:xe, ys:ye)
       case (type_minus)
          a_2d(xs:xe, ys:ye) = b_2d(xs:xe, ys:ye) - c_2d(xs:xe, ys:ye)
       case (type_mult)
          a_2d(xs:xe, ys:ye) = b_2d(xs:xe, ys:ye) * c_2d(xs:xe, ys:ye)
       case (type_divd)
          a_2d(xs:xe, ys:ye) = b_2d(xs:xe, ys:ye) / c_2d(xs:xe, ys:ye)
       end select

       call DMDAVecRestoreArrayF90(A_dm, A, a_2d, ierr)
       call DMDAVecRestoreArrayF90(B_dm, B, b_2d, ierr)
       call DMDAVecRestoreArrayF90(C_dm, C, c_2d, ierr)
    case (3)
       call DMDAVecGetArrayF90(A_dm, A, a_3d, ierr)
       call DMDAVecGetArrayF90(B_dm, B, b_3d, ierr)
       call DMDAVecGetArrayF90(C_dm, C, c_3d, ierr)

       select case(op_type)
       case (type_plus)
          a_3d(xs:xe, ys:ye, zs:ze) = b_3d(xs:xe, ys:ye, zs:ze) + &
               c_3d(xs:xe, ys:ye, zs:ze)
       case (type_minus)
          a_3d(xs:xe, ys:ye, zs:ze) = b_3d(xs:xe, ys:ye, zs:ze) - &
               c_3d(xs:xe, ys:ye, zs:ze)
       case (type_mult)
          a_3d(xs:xe, ys:ye, zs:ze) = b_3d(xs:xe, ys:ye, zs:ze) * &
               c_3d(xs:xe, ys:ye, zs:ze)
       case (type_divd)
          a_3d(xs:xe, ys:ye, zs:ze) = b_3d(xs:xe, ys:ye, zs:ze) / &
               c_3d(xs:xe, ys:ye, zs:ze)
       end select

       call DMDAVecRestoreArrayF90(A_dm, A, a_3d, ierr)
       call DMDAVecRestoreArrayF90(B_dm, B, b_3d, ierr)
       call DMDAVecRestoreArrayF90(C_dm, C, c_3d, ierr)

    end select

    ! write(*, "(A, 100g15.5)") "a_x1=", a_x1(xs:xe,ys:ye,ze:ze)    
    
    call PetscLogEventEnd(ievent,ierr)
  end subroutine 
end module 
