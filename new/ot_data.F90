#:include "type_def.fypp"
#include "type.h"

#define I1 :
#define I2 :,:
#define I3 :,:,:

#define D1 xs:xe
#define D2 xs:xe,ys:ye
#define D3 xs:xe,ys:ye,zs:ze

module ot_data
  use ot_common
  use ot_petsc
  use ot_geom

  implicit none
  
  interface get_local_ptr
     module procedure get_local_ptr1d
     module procedure get_local_ptr2d
     module procedure get_local_ptr3d
  end interface get_local_ptr
  
contains

#include "data/data_rand.F90"
#include "data/data_consts.F90"
#include "data/data_seqs.F90"
  
  subroutine init_data(ierr)
    implicit none
#include "petsc.h"
    integer, intent(out) :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  end subroutine 

  subroutine finalize_data(ierr)
    implicit none
#include "petsc.h"    
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

    call assert(A /= 0, __FILE__, __LINE__, &
         "${fun_name}$ : result vec not allocated!")

    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    call VecGetDM(A, A_dm, ierr)
    call VecGetDM(B(1), B_dm, ierr)

    call petsc_get_corners(A, box)

    xs = box%rx%lower
    ys = box%ry%lower
    zs = box%rz%lower

    xe = box%rx%upper
    ye = box%ry%upper
    ze = box%rz%upper

    call DMDAVecGetArrayF90(A_dm, A, a_3d, ierr)

    a_3d(D3) = 0.

    do i = 1, size(B)
       call DMDAVecGetArrayF90(B_dm, B(i), b_3d, ierr)

       a_3d(D3) = a_3d(D3) ${op[3]}$ &
            (ops_alpha(i) + ops_beta(i) * b_3d(D3))

       call DMDAVecRestoreArrayF90(B_dm, B(i), b_3d, ierr)
    enddo

    !a_3d(D3) = alpha + beta * a_3d(D3)

    call DMDAVecRestoreArrayF90(A_dm, A, a_3d, ierr)

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

    call assert(A /= 0, __FILE__, __LINE__, &
         "${fun_name}$ : result vec not allocated!")
    
    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    call VecGetDM(A, A_dm, ierr)
    call VecGetDM(B, B_dm, ierr)
    
    call petsc_get_corners(A, box)

    xs = box%rx%lower
    ys = box%ry%lower
    zs = box%rz%lower

    xe = box%rx%upper
    ye = box%ry%upper
    ze = box%rz%upper

    call DMDAVecGetArrayF90(A_dm, A, a_3d, ierr)
    call DMDAVecGetArrayF90(B_dm, B(1), b1_3d, ierr)
    call DMDAVecGetArrayF90(B_dm, B(2), b2_3d, ierr)       

    a_3d(D3) =  &
         !alpha + beta * &
         (merge(1.d0, 0.d0, &
         (ops_alpha(1) + ops_beta(1) * b1_3d(D3)) &
         ${op[3]}$ &
         (ops_alpha(2) + ops_beta(2) * b2_3d(D3))))

    call DMDAVecRestoreArrayF90(B_dm, B(1), b1_3d, ierr)
    call DMDAVecRestoreArrayF90(B_dm, B(2), b2_3d, ierr)          
    call DMDAVecRestoreArrayF90(A_dm, A, a_3d, ierr)

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
    
    call assert(A /= 0, __FILE__, __LINE__, &
         "${fun_name}$ : result vec not allocated!")
    
    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    call VecGetDM(A, A_dm, ierr)
    call VecGetDM(B, B_dm, ierr)
    
    !call dm_get_corners(A_dm, xs, ys, zs, xl, yl, zl)
    call petsc_get_corners(A, box)

    xs = box%rx%lower
    ys = box%ry%lower
    zs = box%rz%lower

    xe = box%rx%upper
    ye = box%ry%upper
    ze = box%rz%upper
    
    call DMDAVecGetArrayF90(A_dm, A, a_3d, ierr)
    call DMDAVecGetArrayF90(B_dm, B(1), b1_3d, ierr)

#:if op[0] == 'type_pow'
    a_3d(D3) = &
         (ops_alpha(1) + ops_beta(1) * b1_3d (D3))**args(1)
#:else
    a_3d(D3) = &
         ${op[3]}$(ops_alpha(1) + ops_beta(1) * b1_3d (D3))
#:endif

    call DMDAVecRestoreArrayF90(B_dm, B(1), b1_3d, ierr)
    call DMDAVecRestoreArrayF90(A_dm, A, a_3d, ierr)

    ! print*, "ops_alpha=", ops_alpha(1), "ops_beta=", ops_beta(1)
    ! print*, b1_3d (D3)
    ! call VecView(A, PETSC_VIEWER_STDOUT_WORLD,ierr)
    ! call VecView(B(1), PETSC_VIEWER_STDOUT_WORLD,ierr)       
    
    call PetscLogEventEnd(ievent, ierr)
  end subroutine
  
#:endif
#:endfor

subroutine data_get_sub(dst, src, ref)
  implicit none
#include "petsc.h"
  Vec, intent(out) :: dst   
  Vec, intent(in)  :: src
  DM :: src_dm, dst_dm
  type(ref_info) :: ref
  type(box_info) :: dst_box
  type(box_info) :: dst_local, src_local, dst_local_dst, src_box
  integer :: dst_shape(3), dst_local_shape(3)
  integer :: m, n, k, bm, bn, bk
  integer :: xs, ys, zs, xe, ye, ze
  integer :: xs1, ys1, zs1, xe1, ye1, ze1
  PetscScalar, pointer :: src_data(:,:,:), dst_data(:,:,:)
  integer :: v_shape(3)

  call petsc_get_corners(src, src_local)
  call petsc_get_shape(src_box, src)

  ! if(.not. (dst_box .in. src_local)) then
  !    call disp_box(dst_box, 'dst_box = ')
  !    call disp_box(src_box, 'src_box = ')
  ! end if

  call assert(dst_box .in. src_box, &
       __FILE__, __LINE__, &
       "dst box must be smaller than source box")

  if(ref%ref_index_type_x == 0 .and. &
       ref%ref_index_type_y == 0 .and. ref%ref_index_type_z == 0) then

     call range_to_box(dst_box, ref%range_x, ref%range_y, ref%range_z)

     !dst_local : the coordinates refer to the source data
     dst_local = (src_local .and. dst_box)

     ! call disp_box(src_local, 'src_local = ')
     ! call disp_box(dst_box,   'dst_box = ')
     ! call disp_box(dst_local, 'dst_local = ')

     dst_local_shape = shape(dst_local)
     dst_shape = shape(dst_box)

     call petsc_get_dm(src_dm, src)
     call petsc_slice_dm(dst_dm, src_dm, dst_box)
     call petsc_new3d(dst, dst_dm)

     !call disp(src_local, 'src_local = ')

     call petsc_get_corners(dst, dst_local_dst)

     xs1 = dst_local_dst%rx%lower
     ys1 = dst_local_dst%ry%lower
     zs1 = dst_local_dst%rz%lower
     xe1 = dst_local_dst%rx%upper
     ye1 = dst_local_dst%ry%upper
     ze1 = dst_local_dst%rz%upper

     xs = dst_local%rx%lower
     ys = dst_local%ry%lower
     zs = dst_local%rz%lower
     xe = dst_local%rx%upper
     ye = dst_local%ry%upper
     ze = dst_local%rz%upper

     call get_local_arr(dst, dst_data)
     call get_local_arr(src, src_data)

     dst_data(xs1:xe1,ys1:ye1,zs1:ze1) = &    
          src_data(xs:xe,ys:ye,zs:ze)

     call restore_local_arr(dst, dst_data)
  end if
end subroutine

#:set itype = [['range', 'type(range)'], ['iarr', 'integer, dimension(:)']]
#:for tx in itype
#:for ty in itype
#:for tz in itype
  subroutine data_get_sub_${tx[0]}$_${ty[0]}$_${tz[0]}$ &
       (dst, src, ix, iy, iz)
    implicit none
#include "petsc.h"
    Vec, intent(out) :: dst   
    Vec, intent(in)  :: src
    ${tx[1]}$ :: ix
    ${ty[1]}$ :: iy
    ${tz[1]}$ :: iz    
    DM :: src_dm, dst_dm
    type(box_info) :: dst_local, src_local, dst_local_dst, src_box
    integer :: dst_shape(3), dst_local_shape(3)
    integer :: m, n, k, bm, bn, bk
    integer :: xs, ys, zs, xe, ye, ze
    integer :: xs1, ys1, zs1, xe1, ye1, ze1
    PetscScalar, pointer :: src_data(:,:,:), dst_data(:,:,:)
    integer :: v_shape(3)
    
    ! call petsc_get_corners(src, src_local)
    ! call petsc_get_shape(src_box, src)

    ! ! if(.not. (dst_box .in. src_local)) then
    ! !    call disp_box(dst_box, 'dst_box = ')
    ! !    call disp_box(src_box, 'src_box = ')
    ! ! end if
    
    ! !dst_local : the coordinates refer to the source data
    ! dst_local = (src_local .and. dst_box)
    
    ! ! call disp_box(src_local, 'src_local = ')
    ! ! call disp_box(dst_box,   'dst_box = ')
    ! ! call disp_box(dst_local, 'dst_local = ')
    
    ! dst_local_shape = shape(dst_local)
    ! dst_shape = shape(dst_box)

    ! call petsc_get_dm(src_dm, src)
    ! call petsc_slice_dm(dst_dm, src_dm, dst_box)
    ! call petsc_new3d(dst, dst_dm)
    
    ! !call disp(src_local, 'src_local = ')

    ! call petsc_get_corners(dst, dst_local_dst)
    
    ! xs1 = dst_local_dst%rx%lower
    ! ys1 = dst_local_dst%ry%lower
    ! zs1 = dst_local_dst%rz%lower
    ! xe1 = dst_local_dst%rx%upper
    ! ye1 = dst_local_dst%ry%upper
    ! ze1 = dst_local_dst%rz%upper

    ! xs = dst_local%rx%lower
    ! ys = dst_local%ry%lower
    ! zs = dst_local%rz%lower
    ! xe = dst_local%rx%upper
    ! ye = dst_local%ry%upper
    ! ze = dst_local%rz%upper

    ! call get_local_arr(dst, dst_data)
    ! call get_local_arr(src, src_data)

    ! dst_data(xs1:xe1,ys1:ye1,zs1:ze1) = &    
    !      src_data(xs:xe,ys:ye,zs:ze)
    
    ! call restore_local_arr(dst, dst_data)

  end subroutine
#:endfor
#:endfor
#:endfor
  
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

  subroutine data_set_scalar(data, val, ref)
    implicit none
    Vec :: data
    real(8) :: val
    type(ref_info) :: ref
    PetscScalar, pointer :: local_arr(:,:,:)
    type(box_info) :: local_box, global_box
    type(box_info) :: set_local_box, set_box
    integer :: xs, xe, ys, ye, zs, ze
    integer :: ierr

    if(ref%ref_index_type_x == 0 .and. &
         ref%ref_index_type_y == 0 .and. &
         ref%ref_index_type_z == 0) then

       if(ref%range_x == range_all .and. &
            ref%range_y == range_all .and. &
            ref%range_z == range_all) then
          call VecSet(data, val, ierr)
       else
          call range_to_box(set_box, ref%range_x, ref%range_y, ref%range_z)
          call petsc_get_shape(global_box, local_box, data)
          set_local_box = (local_box .and. set_box)

          call get_local_arr(data, local_arr)
          xs = set_local_box%rx%lower
          ys = set_local_box%ry%lower
          zs = set_local_box%rz%lower
          xe = set_local_box%rx%upper
          ye = set_local_box%ry%upper
          ze = set_local_box%rz%upper
          local_arr(xs:xe,ys:ye,zs:ze) = val
          call restore_local_arr(data, local_arr)
       endif
    endif
  end subroutine

  subroutine data_set_ref_ref(dst, src, dst_ref, src_ref)
    implicit none
#include "petsc.h"
    Vec, intent(inout) :: dst
    Vec, intent(in) :: src
    DM :: src_dm, dst_dm
    type(ref_info), intent(in) :: dst_ref, src_ref
    type(dist_info) :: dst_dist, src_dist, dst_new_dist, src_new_dist
    integer :: ierr
    integer :: xs, ys, zs, xe, ye, ze
    integer :: xs1, ys1, zs1, xe1, ye1, ze1    
    PetscScalar, pointer :: src_ptr(:,:,:)
    PetscScalar, pointer :: src_ptr1d(:)
    PetscScalar, pointer :: dst_ptr(:,:,:)
    type(box_info) :: src_local_box, dst_local_box
    type(box_info) :: src_ref_box,   dst_ref_box
    type(box_info) :: dst_local_ref_box, src_local_ref_box
    integer :: rank
    integer, allocatable :: idx(:)
    integer :: src_shape(3), ref_local_shape(3)
    Vec :: vec_out
    IS :: is_out
    DM :: vec_dm
    
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    call petsc_get_dm(src_dm, src)
    call petsc_get_dm(dst_dm, dst)
    call petsc_get_dist(dst_dist, dst_dm)
    call petsc_get_dist(src_dist, src_dm)

    call petsc_new_dist(dst_new_dist, dst_dist, dst_ref)
    call petsc_new_dist(src_new_dist, src_dist, src_ref)

    ! write(*, "(A, 30I7.0)") "dst_dist%lx=", dst_dist%lx
    ! write(*, "(A, 30I7.0)") "dst_dist%ly=", dst_dist%ly
    ! write(*, "(A, 30I7.0)") "dst_dist%lz=", dst_dist%lz    
    ! write(*, "(A, 30I7.0)") "src_dist%lx=", src_dist%lx
    ! write(*, "(A, 30I7.0)") "src_dist%ly=", src_dist%ly
    ! write(*, "(A, 30I7.0)") "src_dist%lz=", src_dist%lz    

    call petsc_get_corners(dst, dst_local_box)
    call petsc_get_corners(src, src_local_box)

    call range_to_box(src_ref_box, src_ref%range_x, &
         src_ref%range_y, src_ref%range_z)

    call range_to_box(dst_ref_box, dst_ref%range_x, &
         dst_ref%range_y, dst_ref%range_z)

    dst_local_ref_box = (dst_local_box .and. dst_ref_box)
    ! src_local_ref_box = (src_local_box .and. src_ref_box)

    call get_ref_box(src_local_ref_box, src_ref_box, &
         dst_local_ref_box, dst_ref_box)

    ! call disp(src_local_ref_box, 'src_local_ref_box = xxx')    
    ! call disp(src_ref_box, 'src_ref_box = xxx')
    ! call disp(dst_local_ref_box, 'dst_local_ref_box = xxx')
    ! call disp(dst_ref_box, 'dst_ref_box = xxx')    
    
    xs = dst_local_ref_box%rx%lower
    ys = dst_local_ref_box%ry%lower
    zs = dst_local_ref_box%rz%lower
    xe = dst_local_ref_box%rx%upper
    ye = dst_local_ref_box%ry%upper
    ze = dst_local_ref_box%rz%upper       
    
    !only local copy !!
    if(dst_new_dist == src_new_dist) then
       if(is_valid(dst_local_ref_box)) then
          ! call mpi_order_start(MPI_COMM_WORLD,ierr)
          ! write(*, "(A, I3.0)") "rank = ", get_rank(MPI_COMM_WORLD)
          ! call disp(dst_local_box, 'dst_local_box = ')
          ! call disp(src_local_box, 'src_local_box = ')       
          ! call disp(dst_ref_box, 'dst_ref_box = ')
          ! call disp(src_ref_box, 'src_ref_box = ')       
          ! call disp(dst_local_ref_box, 'dst_local_ref_box = ')
          ! call disp(src_local_ref_box, 'src_local_ref_box = ')
          ! call mpi_order_end(MPI_COMM_WORLD, ierr)              

          xs1 = src_local_ref_box%rx%lower
          ys1 = src_local_ref_box%ry%lower
          zs1 = src_local_ref_box%rz%lower
          xe1 = src_local_ref_box%rx%upper
          ye1 = src_local_ref_box%ry%upper
          ze1 = src_local_ref_box%rz%upper       

          call get_local_arr(src, src_ptr)
          call get_local_arr(dst, dst_ptr)

          dst_ptr(xs:xe, ys:ye, zs:ze) = &
               src_ptr(xs1:xe1, ys1:ye1, zs1:ze1)

          call restore_local_arr(dst, dst_ptr)
          call restore_local_arr(src, src_ptr)
       end if
    else
       call petsc_get_shape(src_shape, src)
       
       call box_to_indices(idx, src_local_ref_box, src_shape)

       call ISCreateGeneral(MPI_COMM_WORLD, size(idx), &
            idx, PETSC_COPY_VALUES, is_out, ierr)
       
       call VecGetSubVector(src, is_out, vec_out, ierr)

       !NODE: the index of returned array pointer starts from 1
       call VecGetArrayF90(vec_out, src_ptr1d, ierr)
       
       !call VecView(vec_out, PETSC_VIEWER_STDOUT_WORLD, ierr)
       !call VecGetDM(vec_out, vec_dm, ierr)

       if(is_valid(dst_local_ref_box)) then
          call get_local_arr(dst, dst_ptr)

          ! call disp(dst_local_ref_box, 'dst_local_ref_box=')
          ref_local_shape = shape(dst_local_ref_box)

          ! print*, "ref_local_shape = ", ref_local_shape
          ! print*, "size(src_ptr1d) = ", size(src_ptr1d)
          ! print*, "src_ptr1d = ", reshape(src_ptr1d, ref_local_shape)

          dst_ptr(xs:xe, ys:ye, zs:ze) = &
               reshape(src_ptr1d, ref_local_shape)

          call VecRestoreArrayF90(vec_out, src_ptr1d, ierr)

          call restore_local_arr(dst, dst_ptr)
       end if
       !call get_local_arr(vec_out, src_ptr1d)
       !ref_local_shape = shape(dst_local_ref_box)
       !call restore_local_arr(vec_out, src_ptr1d)
       
       !call VecView(vec_out, PETSC_VIEWER_STDOUT_WORLD, ierr)
       call ISDestroy(is_out, ierr)
       call VecDestroy(vec_out, ierr)
    end if

  end subroutine
end module

#undef D1
#undef D2
#undef D3


