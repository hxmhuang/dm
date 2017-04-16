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
  use ot_type
  use ot_buffer
  use ot_tensor
  
  interface get_local_ptr
     module procedure get_local_ptr1d
     module procedure get_local_ptr2d
     module procedure get_local_ptr3d
  end interface get_local_ptr
  
contains

#include "data/data_rand.F90"
#include "data/data_consts.F90"
#include "data/data_seqs.F90"
  
  
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
  subroutine ${fun_name}$ (A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num, ierr) 
    implicit none
#include "petsc.h"
    type(tensor), intent(inout)  :: A
    type(tensor_ptr), intent(in) :: B(2)
    PetscScalar, intent(in) :: ops_alpha(2), ops_beta(2)
    PetscScalar, intent(in) :: ops_args(args_num)
    integer, intent(in) :: ops_num
    integer, intent(in) :: args_num
    DM :: A_dm, B_dm
    PetscErrorCode, intent(out) ::ierr
    PetscInt :: xs, xl, ys, yl, zs, zl, xe, ye, ze
    PetscScalar, pointer :: a_3d(:,:,:), b1_3d(:,:,:), b2_3d(:,:,:)
    integer :: i
    type(box_info) :: box
    PetscLogEvent :: ievent
    
    call PetscLogEventRegister("${fun_name}$", 0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    print*, "calling function ", "${fun_name}$"

    ! call disp_info(B(1)%ptr, 'B(1)%ptr = ')
    ! call disp_info(B(2)%ptr, 'B(2)%ptr = ')
    
    ! write(*, "(A, Z16.16)"), "B(1)%ptr = ", B(1)%ptr%data
    ! write(*, "(A, Z16.16)"), "B(2)%ptr = ", B(2)%ptr%data
    
    ! if(size(B) < 2) then
    !    print*, &
    !         "Error! the number of operands for ${fun_name}$ &
    !         & is incorrect. size(operands) = ", size(B)
    !    call abort()
    ! end if

    ! call assert(A /= 0, __FILE__, __LINE__, &
    !      "${fun_name}$ : result vec not allocated!")

    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    call petsc_get_corners(A%data, box)

    xs = box%rx%lower;  xe = box%rx%upper
    ys = box%ry%lower;  ye = box%ry%upper
    zs = box%rz%lower;  ze = box%rz%upper

    call get_local_arr(A%data, a_3d)    
    call get_local_arr(B(1)%ptr%data, b1_3d)
    call get_local_arr(B(2)%ptr%data, b2_3d)
    
    a_3d(D3) = &
         (ops_alpha(1) + ops_beta(1) * b1_3d(D3)) &
         ${op[3]}$ &
         (ops_alpha(2) + ops_beta(2) * b2_3d(D3))

    call restore_local_arr(A%data, a_3d)    
    call restore_local_arr(B(1)%ptr%data, b1_3d)
    call restore_local_arr(B(2)%ptr%data, b2_3d)

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
#:set fun_name = 'data_{0}'.format(op[2])
  !> array ${op[2]}$
  subroutine ${fun_name}$ (A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num, ierr) 
    implicit none
#include "petsc.h"
    type(tensor), intent(inout)  :: A
    type(tensor_ptr), intent(in) :: B(2)
    PetscScalar :: ops_alpha(2), ops_beta(2) !, alpha, beta
    !this argument won't be used    
    !type(buffer_r8), intent(in), pointer :: ops_args(args_num)
    PetscScalar, intent(in) :: ops_args(args_num)
    integer, intent(in) :: ops_num
    integer, intent(in) :: args_num    
    PetscErrorCode, intent(out) ::ierr    
    DM :: A_dm, B_dm
    PetscInt :: xs, ys, zs, xe, ye, ze
    PetscScalar, pointer :: a_3d(:,:,:)
    PetscScalar, pointer :: b1_3d(:,:,:), b2_3d(:,:,:)
    type(box_info) :: box
    PetscLogEvent  ::  ievent
    
    call PetscLogEventRegister("${fun_name}$", 0, ievent, ierr)
    call PetscLogEventBegin(ievent, ierr)

    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    call petsc_get_corners(A%data, box)

    xs = box%rx%lower; xe = box%rx%upper
    ys = box%ry%lower; ye = box%ry%upper
    zs = box%rz%lower; ze = box%rz%upper

    call get_local_arr(A%data,        a_3d)
    call get_local_arr(B(1)%ptr%data, b1_3d)
    call get_local_arr(B(2)%ptr%data, b2_3d)

    if(is_scalar(B(1)%ptr)) then
       a_3d(D3) = (merge(1.d0, 0.d0, &
            b1_3d(1,1,1) &
            ${op[3]}$ &
            (ops_alpha(2) + ops_beta(2) * b2_3d(D3))))
    else if(is_scalar(B(2)%ptr)) then
       a_3d(D3) = (merge(1.d0, 0.d0, &
            (ops_alpha(1) + ops_beta(1) * b1_3d(D3)) &
            ${op[3]}$ &
            b2_3d(1,1,1)))
    else
       a_3d(D3) = (merge(1.d0, 0.d0, &
            (ops_alpha(1) + ops_beta(1) * b1_3d(D3)) &
            ${op[3]}$ &
            (ops_alpha(2) + ops_beta(2) * b2_3d(D3))))
    end if
 
    call restore_local_arr(A%data, a_3d)
    call restore_local_arr(B(1)%ptr%data, b1_3d)
    call restore_local_arr(B(2)%ptr%data, b2_3d)
    
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
  subroutine ${fun_name}$ (A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num, ierr) 
    implicit none
#include "petsc.h"
    type(tensor), intent(inout) :: A
    type(tensor_ptr), intent(in)  :: B(ops_num)
    PetscScalar  :: ops_alpha(ops_num), ops_beta(ops_num) !, alpha, beta
    !type(buffer_r8), intent(in), pointer :: ops_args(args_num)
    PetscScalar, intent(in) :: ops_args(args_num)
    integer, intent(in) :: ops_num
    integer, intent(in) :: args_num
    DM :: A_dm, B_dm
    PetscErrorCode, intent(out) ::ierr
    PetscLogEvent  ::  ievent
    PetscInt :: xs, xl, ys, yl, zs, zl, xe, ye, ze
    PetscScalar, pointer :: a_3d(:,:,:), b1_3d(:,:,:)
    integer :: i
    PetscScalar :: power
    type(box_info) :: box

    call PetscLogEventRegister("${fun_name}$", 0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    if(size(B) /= 1) then
       print*, &
            "Error! the number of operands for ${fun_name}$ &
            & is incorrect. size(operands) = ", size(B)
       call abort()
    end if
    
    ! call assert(A /= 0, __FILE__, __LINE__, &
    !      "${fun_name}$ : result vec not allocated!")

    ! we have to make sure the dm of A, B, C are exactly same.
    ! or we can not just use the local data
    ! call check_dm()
    
    !call dm_get_corners(A_dm, xs, ys, zs, xl, yl, zl)
    call petsc_get_corners(A%data, box)

    xs = box%rx%lower; xe = box%rx%upper
    ys = box%ry%lower; ye = box%ry%upper
    zs = box%rz%lower; ze = box%rz%upper
    
    call get_local_arr(A%data, a_3d)
    call get_local_arr(B(1)%ptr%data, b1_3d)
    
#:if op[0] == 'type_pow'
    power = ops_args(1)
    a_3d(D3) = &
         (ops_alpha(1) + ops_beta(1) * b1_3d (D3))**power
#:else
    a_3d(D3) = &
         ${op[3]}$(ops_alpha(1) + ops_beta(1) * b1_3d (D3))
#:endif

    call restore_local_arr(A%data, a_3d)    
    call restore_local_arr(B(1)%ptr%data, b1_3d)
    
    ! print*, "ops_alpha=", ops_alpha(1), "ops_beta=", ops_beta(1)
    ! print*, b1_3d (D3)
    ! call VecView(A, PETSC_VIEWER_STDOUT_WORLD,ierr)
    ! call VecView(B(1), PETSC_VIEWER_STDOUT_WORLD,ierr)       
    
    call PetscLogEventEnd(ievent, ierr)
  end subroutine
  
#:endif
#:endfor

  subroutine data_matmul(C, A, B)
    implicit none
    Vec, intent(in) :: A, B
    Vec, intent(out) :: C
    
  end subroutine
  
subroutine slice_tensor(dst, src, ref)
  implicit none
#include "petsc.h"
  type(tensor), intent(inout) :: dst   
  type(tensor), intent(in)  :: src
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
  Vec :: vec_dst
  
  call petsc_get_corners(src%data, src_local)
  call petsc_get_shape(src_box, src%data)

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

     call petsc_get_dm(src_dm, src%data)
     call petsc_slice_dm(dst_dm, src_dm, dst_box)
     call petsc_create3d(vec_dst, dst_dm)
     call tensor_bind_data(dst, vec_dst)
     
     call petsc_destroy_dm(dst_dm)
     
     !call disp(src_local, 'src_local = ')

     call petsc_get_corners(dst%data, dst_local_dst)

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

     call get_local_arr(dst%data, dst_data)
     call get_local_arr(src%data, src_data)

     dst_data(xs1:xe1,ys1:ye1,zs1:ze1) = &    
          src_data(xs:xe,ys:ye,zs:ze)

     call restore_local_arr(dst%data, dst_data)
     call restore_local_arr(src%data, src_data)
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
  
 
  !> set a scalar to tensor in a given range
  subroutine data_set_scalar(data, val, ref)
    implicit none
    type(tensor), intent(in) :: data
    real(8) :: val
    type(ref_info) :: ref
    PetscScalar, pointer :: local_arr(:,:,:)
    type(box_info) :: local_box, global_box
    type(box_info) :: set_local_box, set_box
    integer :: xs, xe, ys, ye, zs, ze
    integer :: ierr
    PetscLogEvent  ::  ievent

    call PetscLogEventRegister("data_set_scalar", 0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    ! call mpi_order_start(MPI_COMM_WORLD, ierr)
    ! call disp(ref, 'ref = ')
    ! call mpi_order_end(MPI_COMM_WORLD, ierr)
    
    if(ref%ref_index_type_x == 0 .and. &
         ref%ref_index_type_y == 0 .and. &
         ref%ref_index_type_z == 0) then

       if(ref%range_x   == range_all .and. &
            ref%range_y == range_all .and. &
            ref%range_z == range_all) then
          call VecSet(data%data, val, ierr)
       else
          call range_to_box(set_box, ref%range_x, ref%range_y, ref%range_z)
          call petsc_get_shape(global_box, local_box, data%data)
          set_local_box = (local_box .and. set_box)

          ! call mpi_order_start(MPI_COMM_WORLD, ierr)
          ! call disp(set_local_box, 'set_local_box = ')
          ! call mpi_order_end(MPI_COMM_WORLD, ierr)
          
          call get_local_arr(data%data, local_arr)
          xs = set_local_box%rx%lower
          ys = set_local_box%ry%lower
          zs = set_local_box%rz%lower
          xe = set_local_box%rx%upper
          ye = set_local_box%ry%upper
          ze = set_local_box%rz%upper
          local_arr(xs:xe,ys:ye,zs:ze) = val
          call restore_local_arr(data%data, local_arr)
       endif
    endif

    call PetscLogEventBegin(ievent,ierr)
  end subroutine

  subroutine data_set_ref_ref(dst, src, dst_ref, src_ref)
    implicit none
#include "petsc.h"
    type(tensor), intent(inout) :: dst
    type(tensor), intent(in) :: src
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
    integer, allocatable :: idx(:,:,:)
    integer :: src_shape(3), ref_local_shape(3)
    Vec :: vec_out
    IS :: is_out
    DM :: vec_dm
    integer, allocatable :: vec_dist(:)
    integer :: i
    
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    call petsc_get_dm(src_dm, src%data)
    call petsc_get_dm(dst_dm, dst%data)
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

    call petsc_get_corners(dst%data, dst_local_box)
    call petsc_get_corners(src%data, src_local_box)

    call range_to_box(src_ref_box, src_ref%range_x, &
         src_ref%range_y, src_ref%range_z)

    call range_to_box(dst_ref_box, dst_ref%range_x, &
         dst_ref%range_y, dst_ref%range_z)

    dst_local_ref_box = (dst_local_box .and. dst_ref_box)
    ! src_local_ref_box = (src_local_box .and. src_ref_box)

    call get_ref_box(src_local_ref_box, src_ref_box, &
         dst_local_ref_box, dst_ref_box)

    ! call mpi_order_start(MPI_COMM_WORLD, ierr)
    ! print*, "--------------------------------"
    ! call disp(src_ref_box, 'src_ref_box = xxx')
    ! call disp(dst_ref_box, 'dst_ref_box = xxx')        
    ! call disp(src_local_ref_box, 'src_local_ref_box = xxx')
    ! call disp(dst_local_ref_box, 'dst_local_ref_box = xxx')    
    ! call mpi_order_end(MPI_COMM_WORLD, ierr)
    
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

          call get_local_arr(src%data, src_ptr)
          call get_local_arr(dst%data, dst_ptr)

          dst_ptr(xs:xe, ys:ye, zs:ze) = &
               src_ptr(xs1:xe1, ys1:ye1, zs1:ze1)

          call restore_local_arr(dst%data, dst_ptr)
          call restore_local_arr(src%data, src_ptr)
       end if
    else
       call petsc_get_shape(src_shape, src%data)

       ! call VecView(src%data, PETSC_VIEWER_STDOUT_WORLD, ierr)
       
       !call box_to_indices1(idx, src_local_ref_box, src_shape)

       call petsc_get_vec_dist(vec_dist, src%data)
       call box_to_indices2(idx, src_dist, src_local_ref_box, vec_dist)
       
       ! call mpi_order_start(MPI_COMM_WORLD, ierr)       
       ! print*, "idx = ", idx
       ! print*, 'src_shape = ', src_shape
       ! call disp(src_local_ref_box, 'src_local_ref_box = ')
       ! call mpi_order_end(MPI_COMM_WORLD, ierr)

       ! do i = 1, size(idx)
       !    idx(i) = i - 1
       ! enddo
       
       call ISCreateGeneral(MPI_COMM_WORLD, size(idx), &
            idx, PETSC_COPY_VALUES, is_out, ierr)
       
       call VecGetSubVector(src%data, is_out, vec_out, ierr)

       !call mpi_order_start(MPI_COMM_WORLD, ierr)
       !print*, "gshape = ", src_shape
       !call disp(src_local_ref_box, 'box = ')
       
       !print*, "src_ptr1d = ", src_ptr1d
       !call mpi_order_end(MPI_COMM_WORLD, ierr)    
       
       ! call VecView(vec_out, PETSC_VIEWER_STDOUT_WORLD, ierr)
       
       !NODE: the index of returned array pointer starts from 1
       call VecGetArrayF90(vec_out, src_ptr1d, ierr)

       ! if(get_rank() == 0) then
       !    print*, "rank = ", get_rank(), " data = ", src_ptr1d
       !    print*, "idx = ", idx
       !    print*, 'src_shape = ', src_shape
       !    call disp(src_local_ref_box, 'src_local_ref_box = ')
       ! end if
    
       !call VecGetDM(vec_out, vec_dm, ierr)

       if(is_valid(dst_local_ref_box)) then
          call get_local_arr(dst%data, dst_ptr)

          ! call disp(dst_local_ref_box, 'dst_local_ref_box=')
          ref_local_shape = shape(dst_local_ref_box)

          ! print*, "ref_local_shape = ", ref_local_shape
          ! print*, "size(src_ptr1d) = ", size(src_ptr1d)
          ! print*, "src_ptr1d = ", reshape(src_ptr1d, ref_local_shape)

          dst_ptr(xs:xe, ys:ye, zs:ze) = &
               reshape(src_ptr1d, ref_local_shape)

          call VecRestoreArrayF90(vec_out, src_ptr1d, ierr)

          call restore_local_arr(dst%data, dst_ptr)
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


