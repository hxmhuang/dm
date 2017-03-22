
#define I1 :
#define I2 :,:
#define I3 :,:,:

#define D1 xs:xe
#define D2 xs:xe,ys:ye
#define D3 xs:xe,ys:ye,zs:ze

module ot_data
  use ot_common
  use petsc_helper
#:include "type_def.fypp"
  implicit none
#include "node_type.h"
  
  interface get_local_ptr
     module procedure get_local_ptr1d
     module procedure get_local_ptr2d
     module procedure get_local_ptr3d
  end interface get_local_ptr

  
contains

#include "data/data_rand.F90"
#include "data/data_consts.F90"
  
  subroutine init_data(ierr)
    implicit none
#include "petsc.h"
    integer, intent(out) :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  end subroutine 

  subroutine dm_finalize(ierr)
    integer, intent(out) :: ierr
    
    call PetscFinalize(ierr)
    
  end subroutine 

  subroutine data_duplicate(dst, src, ierr)
    implicit none
#include "petsc.h"
    Vec, intent(out) :: dst
    Vec, intent(in)  :: src
    ! DM, intent(out)  :: dst_dm
    ! DM, intent(in)   :: src_dm
    integer, intent(out) :: ierr
    
    !call DMClone(src_dm, dst_dm, ierr)

    ! write(*, "(A, Z16.16)"), "src=", src
    ! write(*, "(A, Z16.16)"), "dst=", dst
    
    !call VecView(src,PETSC_VIEWER_STDOUT_WORLD,ierr)    

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
    call DMDAGetCorners(ada,xs, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, &
         xl,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr)
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

  elemental function rcp(a) result (res)
    implicit none
    real(8), intent(in) :: a
    real(8) :: res
    res = 1.0 / a
  end function

  !basic arithmatic operations
#:for op in L
#:if op[4] == 'A'
#:set fun_name = 'data_' + op[2]
  !> array ${op[2]}$
  subroutine ${fun_name}$ (A, B, ops_alpha, ops_beta, args, ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B(:)
    PetscScalar :: ops_alpha(:), ops_beta(:) !, alpha, beta
    DM :: A_dm, B_dm
    PetscErrorCode, intent(out) ::ierr
    PetscInt :: m, n, k
    integer  :: num_dim
    PetscLogEvent            ::  ievent
    PetscInt :: xs, xl, ys, yl, zs, zl, xe, ye, ze
    PetscScalar, pointer :: a_1d(:),     b_1d(:)
    PetscScalar, pointer :: a_2d(:,:),   b_2d(:,:)
    PetscScalar, pointer :: a_3d(:,:,:), b_3d(:,:,:)
    real(8), intent(in) :: args(10)
    integer :: i
    type(box_info) :: box
    
    if(size(B) < 2) then
       print*, &
            "Error! the number of operands for ${fun_name}$ &
            & is incorrect. size(operands) = ", size(B)
       call abort()
    end if

    call PetscLogEventRegister("${fun_name}$", 0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    if(A == 0) &
         call VecDuplicate(B(1), A, ierr)

    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    call VecGetDM(A, A_dm, ierr)
    call VecGetDM(B(1), B_dm, ierr)

    call dm_get_corners(A_dm, box)

    xs = box%starts(1)
    ys = box%starts(2)
    zs = box%starts(3)

    xe = box%ends(1)
    ye = box%ends(2)
    ze = box%ends(3)

    ! xe = xs + xl - 1
    ! ye = ys + yl - 1
    ! ze = zs + zl - 1
    
    call dm_get_dim(A_dm, num_dim)

    select case (num_dim)
#:for d in [1, 2, 3]
    case (${d}$)
       
       call DMDAVecGetArrayF90(A_dm, A, a_${d}$d, ierr)

       a_${d}$d(D${d}$) = 0.
       
       do i = 1, size(B)
          call DMDAVecGetArrayF90(B_dm, B(i), b_${d}$d, ierr)

          a_${d}$d(D${d}$) = a_${d}$d(D${d}$) ${op[3]}$ &
               (ops_alpha(i) + ops_beta(i) * b_${d}$d(D${d}$))
          
          call DMDAVecRestoreArrayF90(B_dm, B(i), b_${d}$d, ierr)
       enddo

       !a_${d}$d(D${d}$) = alpha + beta * a_${d}$d(D${d}$)
       
       call DMDAVecRestoreArrayF90(A_dm, A, a_${d}$d, ierr)
#:endfor
    end select

    !call VecView(A, PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    call PetscLogEventEnd(ievent, ierr)
  end subroutine
#:endif  
#:endfor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!! Compare subrutines
!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#:for op in L
#:if op[4] == 'B'
#:set fun_name = 'data_' + op[2]  
  !> array ${op[2]}$
  subroutine ${fun_name}$ (A, B, ops_alpha, ops_beta, args, ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in) :: B(:)
    PetscScalar :: ops_alpha(:), ops_beta(:)!, alpha, beta
    DM :: A_dm, B_dm
    PetscErrorCode, intent(out) ::ierr
    PetscInt :: m, n, k
    integer  :: num_dim
    PetscLogEvent            ::  ievent
    PetscInt :: xs, xl, ys, yl, zs, zl, xe, ye, ze
    PetscScalar, pointer :: a_1d(:),     b1_1d(:),     b2_1d(:)    
    PetscScalar, pointer :: a_2d(:,:),   b1_2d(:,:),   b2_2d(:,:)  
    PetscScalar, pointer :: a_3d(:,:,:), b1_3d(:,:,:), b2_3d(:,:,:)
    integer :: i
    real(8), intent(in) :: args(10)
    type(box_info) :: box
    
    if(size(B) /= 2) then
       print*, &
            "Error! the number of operands for ${fun_name}$ &
            & is incorrect. size(operands) = ", size(B)
       call abort()
    end if

    call PetscLogEventRegister("${fun_name}$", 0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    if(A == 0) &
         call VecDuplicate(A, B(1), ierr)
    
    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    call VecGetDM(A, A_dm, ierr)
    call VecGetDM(B, B_dm, ierr)
    
    call dm_get_corners(A_dm, box)

    xs = box%starts(1)
    ys = box%starts(2)
    zs = box%starts(3)

    xe = box%ends(1)
    ye = box%ends(2)
    ze = box%ends(3)

    ! xe = xs + xl - 1
    ! ye = ys + yl - 1
    ! ze = zs + zl - 1
    
    call dm_get_dim(A_dm, num_dim)

    select case (num_dim)
#:for d in [1, 2, 3]
    case (${d}$)
       
       call DMDAVecGetArrayF90(A_dm, A, a_${d}$d, ierr)
       call DMDAVecGetArrayF90(B_dm, B(1), b1_${d}$d, ierr)
       call DMDAVecGetArrayF90(B_dm, B(2), b2_${d}$d, ierr)       

       a_${d}$d(D${d}$) =  &
 !           alpha + beta * &
            (merge(1.d0, 0.d0, &
            (ops_alpha(1) + ops_beta(1) * b1_${d}$d(D${d}$)) &
            ${op[3]}$ &
            (ops_alpha(2) + ops_beta(2) * b2_${d}$d(D${d}$))))
              
          call DMDAVecRestoreArrayF90(B_dm, B(1), b1_${d}$d, ierr)
          call DMDAVecRestoreArrayF90(B_dm, B(2), b2_${d}$d, ierr)          
          call DMDAVecRestoreArrayF90(A_dm, A, a_${d}$d, ierr)
#:endfor
    end select
    
    call PetscLogEventEnd(ievent, ierr)
  end subroutine
#:endif  
#:endfor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!! Unary operations
!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#:for op in L
#:if op[4] == 'C'
#:set fun_name = 'data_' + op[2]
  !> array ${op[2]}$
  subroutine ${fun_name}$ (A, B, ops_alpha, ops_beta, args, ierr) 
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: A
    Vec, intent(in)  :: B(:)
    PetscScalar  :: ops_alpha(:), ops_beta(:) !, alpha, beta
    DM :: A_dm, B_dm
    PetscErrorCode, intent(out) ::ierr
    PetscInt :: m, n, k
    integer  :: num_dim
    PetscLogEvent            ::  ievent
    PetscInt :: xs, xl, ys, yl, zs, zl, xe, ye, ze
    PetscScalar, pointer :: a_1d(:),     b1_1d(:)
    PetscScalar, pointer :: a_2d(:,:),   b1_2d(:,:)
    PetscScalar, pointer :: a_3d(:,:,:), b1_3d(:,:,:)
    integer :: i
    real(8) :: args(10)
    type(box_info) :: box
    
    if(size(B) /= 1) then
       print*, &
            "Error! the number of operands for ${fun_name}$ &
            & is incorrect. size(operands) = ", size(B)
       call abort()
    end if

    call PetscLogEventRegister("${fun_name}$", 0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    if(A == 0) &
         call VecDuplicate(A, B(1), ierr)
    
    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    call VecGetDM(A, A_dm, ierr)
    call VecGetDM(B, B_dm, ierr)
    
    !call dm_get_corners(A_dm, xs, ys, zs, xl, yl, zl)
    call dm_get_corners(A_dm, box)

    xs = box%starts(1)
    ys = box%starts(2)
    zs = box%starts(3)

    xe = box%ends(1)
    ye = box%ends(2)
    ze = box%ends(3)
    
    ! xe = xs + xl - 1
    ! ye = ys + yl - 1
    ! ze = zs + zl - 1
    
    call dm_get_dim(A_dm, num_dim)

    select case (num_dim)
#:for d in [1, 2, 3]
    case (${d}$)
       
       call DMDAVecGetArrayF90(A_dm, A, a_${d}$d, ierr)
       call DMDAVecGetArrayF90(B_dm, B(1), b1_${d}$d, ierr)

#:if op[0] == 'type_pow'
       a_${d}$d(D${d}$) = &
!            alpha + beta * &
            ((ops_alpha(1) + ops_beta(1) * b1_${d}$d (D${d}$))**args(1))
#:else
       a_${d}$d(D${d}$) = &
!            alpha + beta * &
            (${op[3]}$(ops_alpha(1) + ops_beta(1) * b1_${d}$d (D${d}$)))
#:endif
       
       call DMDAVecRestoreArrayF90(B_dm, B(1), b1_${d}$d, ierr)
       call DMDAVecRestoreArrayF90(A_dm, A, a_${d}$d, ierr)

       ! print*, "ops_alpha=", ops_alpha(1), "ops_beta=", ops_beta(1)
       ! print*, b1_${d}$d (D${d}$)
       ! call VecView(A, PETSC_VIEWER_STDOUT_WORLD,ierr)
       ! call VecView(B(1), PETSC_VIEWER_STDOUT_WORLD,ierr)       
#:endfor
    end select
    
    call PetscLogEventEnd(ievent, ierr)
  end subroutine
  
#:endif
#:endfor

  subroutine data_get_arr(src, box, arr)
    implicit none
    Vec :: src
    type(box_info) :: box
    real(8), intent(out), allocatable :: arr(:)


  end subroutine

!   subroutine data_get_sub(dst, src, box)
!     implicit none
!     Vec :: src
!     type(box_info) :: src_box
!     type(box_info) :: box
!     type(box_info) :: new_local
!     Vec, intent(out) :: dst
!     integer :: src_dim
!     integer :: xs, ys, zs, xe, ye, ze
!     PetscScalar, pointer :: src_data1(:), src_data2(:,:), src_data3(:,:,:)
!     PetscScalar, pointer :: dst_data1(:), dst_data2(:,:), dst_data3(:,:,:)    
    
!     call vec_get_corners(src, src_box)
!     call vec_get_dim(src, src_dim)
    
!     dst_local = (src_box .and. box)

!     dst_local_shape = shape(dst_local)
!     dst_shape = shape(box)
    
!     select case(src_dim)
! #:for d in [1,2,3]
!     case(${d}$)
       
! #:if d <= 1
!        xs = dst_local%starts(1)
!        xe = dst_local%ends(1)
! #:elif d <=2
!        ys = dst_local%starts(2)
!        zs = dst_local%ends(2)
! #:else
!        zs = dst_local%starts(3)
!        ze = dst_local%ends(3)
! #:endif
       
!        call get_local_arr(src, src_data${d}$)
!        allocate(dst_data(D${d}))
!        dst_data${d}$(D${d}$) = src_data${d}$(D${d}$)
       
! #:endfor
!     end select
    
!   end subroutine

  ! subroutine data_create(data, box, local_box)
  !   implicit none
    
  ! end subroutine
  
  ! subroutine data_assign_scalar(dst, dst_box, v)
  !   implicit none
  !   Vec, intent(inout) :: dst
  !   type(box_info) :: dst_box
  !   PetscScalar :: v

  ! end subroutine
  
  ! subroutine data_assign_ref_to_ref(dst, dst_box, src, src_box)
  !   implicit none
  !   Vec, intent(inout) :: dst
  !   type(box_info) :: dst_box, src_box
  !   Vec, intent(in) :: src


  ! end subroutine
  
end module

#undef D1
#undef D2
#undef D3


