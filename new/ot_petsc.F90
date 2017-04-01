
#define I1 :
#define I2 :,:
#define I3 :,:,:

module ot_petsc
  use ot_geom
  use ot_print
  use ot_common
  
#:for d in [1,2,3]
  interface get_local_arr
     module procedure get_local_arr${d}$
  end interface get_local_arr
#:endfor

#:for d in [1,2,3]
  interface restore_local_arr
     module procedure restore_local_arr${d}$
  end interface restore_local_arr
#:endfor

  interface petsc_new3d
     module procedure petsc_new3d_by_shape
     module procedure petsc_new3d_by_dm
  end interface

  interface petsc_get_shape
     module procedure petsc_get_shape1
     module procedure petsc_get_shape2
     module procedure petsc_get_shape3
     module procedure petsc_get_shape4     
  end interface

  interface petsc_update_dist
     module procedure petsc_update_dist1
     module procedure petsc_update_dist2     
  end interface petsc_update_dist
  
contains
    subroutine vec_get_dim(A, num_dim)
    implicit none
#include "petsc.h"
    Vec :: A
    DM :: DM_A
    integer :: ierr
    integer, intent(out) :: num_dim
    
    call VecGetDM(A, DM_A, ierr)

    call DMDAGetInfo(DM_A, num_dim, &
         PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,ierr)
  end subroutine

  subroutine dm_get_dim(DM_A, num_dim)
    implicit none
#include "petsc.h"
    DM :: DM_A
    integer :: ierr
    integer, intent(out) :: num_dim
    
    call DMDAGetInfo(DM_A, num_dim, &
         PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER,ierr)
  end subroutine

  subroutine petsc_get_dm(data_dm, data)
    implicit none
#include "petsc.h"
    Vec, intent(in) :: data
    DM, intent(out) :: data_dm
    integer :: ierr

    call VecGetDM(data, data_dm, ierr)
    CHKERRQ(ierr)
  end subroutine
  
  subroutine petsc_get_corners(A, b)
    implicit none
#include "petsc.h"    
    Vec :: A
    DM  :: DM_A
    integer :: xs, ys, zs, xl, yl, zl
    type(box_info), intent(out) :: b
    integer :: ierr
    
    call VecGetDM(A, DM_A, ierr)
    call DMDAGetCorners(DM_A, xs, ys, zs, xl, yl, zl,ierr)

    b%rx = r(xs, xs + xl - 1)
    b%ry = r(ys, ys + yl - 1)
    b%rz = r(zs, zs + zl - 1)    
  end subroutine

  subroutine petsc_get_dist(dist, data_dm)
    implicit none
#include "petsc.h"
    DM, intent(in) :: data_dm
    type(dist_info), intent(out) :: dist
    integer :: ierr
    integer :: nx, ny, nz !number of processor in each dimension
    
    call DMDAGetInfo(data_dm, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         nx, ny, nz, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

    CHKERRQ(ierr)
    
    allocate(dist%lx(nx), dist%ly(ny), dist%lz(nz))
    
    call DMDAGetOwnershipRanges(data_dm, dist%lx, dist%ly, dist%lz, ierr)
    CHKERRQ(ierr)
  end subroutine

  
  subroutine petsc_get_shape1(res, data)
    implicit none
#include "petsc.h"
    Vec :: data
    DM :: data_dm
    integer, intent(out) :: res(3)
    integer :: dx, dy, dz, ierr

    call VecGetDM(data, data_dm, ierr)
    CHKERRQ(ierr)    
    call DMDAGetInfo(data_dm, PETSC_NULL_INTEGER, dx, dy, dz, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr)

    ! call DMDAGetInfo(data_dm, 0, dx, dy, dz, &
    !      0, 0, 0, 0, 0, 0, 0, 0, 0, ierr)
    
    CHKERRQ(ierr)
    res(1) = dx
    res(2) = dy
    res(3) = dz
  end subroutine

  subroutine petsc_get_shape2(res, data)
    implicit none
#include "petsc.h"    
    Vec :: data
    DM :: data_dm
    integer :: dx, dy, dz, ierr
    type(box_info), intent(out) :: res

    call VecGetDM(data, data_dm, ierr)
    CHKERRQ(ierr)    
    call DMDAGetInfo(data_dm, PETSC_NULL_INTEGER, dx, dy, dz, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr)
    CHKERRQ(ierr)
    
    ! res%starts = 0
    ! res%ends = (/dx,dy,dz/) - 1
    res%rx = r(0, dx-1)
    res%ry = r(0, dy-1)
    res%rz = r(0, dz-1)    
  end subroutine

  !> get global shape of the data and the local corners
  subroutine petsc_get_shape3(global_shape, local_box, data)
    implicit none
#include "petsc.h"    
    Vec :: data
    DM  :: data_dm
    integer :: dx, dy, dz, ierr
    integer, intent(out) :: global_shape(3)
    type(box_info), intent(out) :: local_box
    integer :: xs, ys, zs, xl, yl, zl
    
    call VecGetDM(data, data_dm, ierr)

    CHKERRQ(ierr)
    call DMDAGetInfo(data_dm, PETSC_NULL_INTEGER, dx, dy, dz, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr)

    CHKERRQ(ierr)
    global_shape = (/dx, dy, dz/)

    call DMDAGetCorners(data_dm, xs, ys, zs, xl, yl, zl,ierr)
    CHKERRQ(ierr)
    
    local_box%rx = r(xs, xs + xl - 1)
    local_box%ry = r(ys, ys + yl - 1)
    local_box%rz = r(zs, zs + zl - 1)    

  end subroutine

  !> get global shape of the data and the local corners
  subroutine petsc_get_shape4(global_box, local_box, data)
    implicit none
#include "petsc.h"    
    Vec :: data
    DM :: data_dm
    integer :: dx, dy, dz, ierr
    integer :: xs, ys, zs, xl, yl, zl
    type(box_info), intent(out) :: global_box, local_box
    
    call VecGetDM(data, data_dm, ierr)
    CHKERRQ(ierr)
    call DMDAGetInfo(data_dm, PETSC_NULL_INTEGER, dx, dy, dz, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr)

    CHKERRQ(ierr)
    
    global_box%rx = r(0, dx-1)
    global_box%ry = r(0, dy-1)
    global_box%rz = r(0, dz-1)    
    
    call DMDAGetCorners(data_dm, xs, ys, zs, xl, yl, zl,ierr)
    CHKERRQ(ierr)
    
    local_box%rx = r(xs, xs + xl - 1)
    local_box%ry = r(ys, ys + yl - 1)
    local_box%rz = r(zs, zs + zl - 1)    
  end subroutine

#:for d in [1,2,3]
  subroutine get_local_arr${d}$(data, arr)
    implicit none
#include "petsc.h"    
    Vec :: data
    DM :: data_dm
    PetscScalar, pointer, intent(out) :: arr(I${d}$)
    integer :: ierr
    
    call VecGetDM(data, data_dm, ierr)
    call DMDAVecGetArrayF90(data_dm, data, arr, ierr)
  end subroutine

  subroutine restore_local_arr${d}$(data, arr)
    implicit none
#include "petsc.h"
    Vec :: data
    DM  :: data_dm
    PetscScalar, pointer, intent(inout) :: arr(I${d}$)
    integer :: ierr
    
    call VecGetDM(data, data_dm, ierr)
    call DMDAVecRestoreArrayF90(data_dm, data, arr, ierr)
  end subroutine  
#:endfor

  subroutine petsc_new3d_by_shape(arr, m, n, k, bm, bn, bk)
    implicit none
#include "petsc.h"    
    integer :: m, n, k, bm, bn, bk
    Vec, intent(out) :: arr
    DM :: arr_dm
    integer :: ierr, xs, ys, zs, xl, yl, zl
    
    call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,   &
         DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,          &
         DMDA_STENCIL_STAR, m,n,k, bm, bn, bk,       &
         1, 0,                                       &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,      &
         PETSC_NULL_INTEGER, arr_dm, ierr)

    CHKERRQ(ierr)
    
    call DMDAGetCorners(arr_dm,xs,ys,zs, xl,yl,zl,ierr)
    call DMCreateGlobalVector(arr_dm, arr, ierr)
  end subroutine

  subroutine petsc_new3d_by_dm(arr, arr_dm)
    implicit none
#include "petsc.h"    
    integer :: m, n, k, bm, bn, bk
    Vec, intent(out) :: arr
    DM, intent(in) :: arr_dm
    integer :: ierr
    
    call DMCreateGlobalVector(arr_dm, arr, ierr)

    CHKERRQ(ierr)    
  end subroutine
  
  subroutine petsc_data_info(data, dx, dy, dz)
    implicit none
#include "petsc.h"
    Vec, intent(in) :: data
    DM :: data_dm
    integer, intent(out) :: dx, dy, dz
    integer :: ierr
    
    call VecGetDM(data, data_dm, ierr)
    call DMDAGetInfo(data_dm, PETSC_NULL_INTEGER, dx, dy, dz, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr)
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> update distribution with a new range from a to b
  subroutine petsc_update_dist1(ll, a, b)
    implicit none
    integer, intent(inout) :: ll(:)
    integer :: a, b
    integer :: cum, i, p
    
    cum = 0
    do i = 0, size(ll) - 1
       p = min(b, cum + ll(i+1) - 1) - max(a, cum) + 1
       cum = cum + ll(i + 1)
       if(p < 0) p = 0
       ll(i+1) = p
    enddo
  end subroutine

  subroutine petsc_update_dist2(ll, r)
    implicit none
    integer, intent(inout) :: ll(:)
    type(range) :: r
    integer :: cum, i, p, a, b

    a = r%lower
    b = r%upper
    call petsc_update_dist(ll, a, b)
  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> create a new distribution from src_dist
  !> according to the reference info
  subroutine petsc_new_dist(dst_dist, src_dist, ref)
    implicit none
    type(dist_info), intent(out) :: dst_dist
    type(dist_info), intent(in)  :: src_dist
    type(ref_info) :: ref
    integer :: a, b
    integer :: cum, i, p
    integer :: nx, ny, nz

    nx = size(src_dist%lx)
    ny = size(src_dist%ly)
    nz = size(src_dist%lz)    

    allocate(dst_dist%lx(nx), &
         dst_dist%ly(ny), dst_dist%lz(nz))

    dst_dist = src_dist
    
    call assert(ref%ref_index_type_x == 0 .and. &
         ref%ref_index_type_y == 0 .and. &
         ref%ref_index_type_z == 0, __FILE__, __LINE__, &
         "Only support for range index now.")

    call petsc_update_dist(dst_dist%lx, ref%range_x)
    call petsc_update_dist(dst_dist%ly, ref%range_y)
    call petsc_update_dist(dst_dist%lz, ref%range_z)    
  end subroutine
  
  ! subroutine petsc_update_dist(ll, ix, local)
  !   implicit none
  !   integer, intent(inout) :: ll(:)
  !   integer, intent(in) :: ix(:)
  !   integer, allocatable, dimension(:), intent(out) :: local
  !   integer :: cum, i, p
    
  !   cum = 0
  !   do i = 0, size(ll) - 1
  !      p = min(b, cum + ll(i+1) - 1) - max(a, cum) + 1
  !      cum = cum + ll(i + 1)
  !      if(p < 0) p = 0
  !      ll(i+1) = p
  !   enddo
  ! end subroutine

  
  subroutine petsc_slice_dm(dst, src, box)
    implicit none
#include  "petsc.h"
    
    DM, intent(out) :: dst
    DM, intent(in) :: src
    type(box_info) :: box
    integer :: dim, gx, gy, gz, nx, ny, nz, i
    integer, allocatable :: lx(:), ly(:), lz(:)
    integer :: cum, ierr, p, xs, xe, ys, ye, zs, ze
    character(len=512) :: buffer

    call DMDAGetInfo(src, dim, gx, gy, gz, &
         nx, ny, nz, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

    allocate(lx(nx), ly(ny), lz(nz))
    
    call DMDAGetOwnershipRanges(src, lx, ly, lz, ierr)
    !print*,"lx=", lx, "ly=",ly, "lz=", lz
    !call PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT, ierr)
    
    xs = box%rx%lower
    ys = box%ry%lower
    zs = box%rz%lower

    xe = box%rx%upper
    ye = box%ry%upper
    ze = box%rz%upper
    
    ! cum = 0
    ! do i = 0, nx - 1
    !    p = min(xe, cum + lx(i+1) - 1) - max(xs, cum) + 1
    !    cum = cum + lx(i + 1)
    !    if(p < 0) p = 0
    !    lx(i+1) = p
    ! enddo
    
    ! cum = 0
    ! do i = 0, ny - 1
    !    p = min(ye, cum + ly(i+1) - 1) - max(ys, cum) + 1
    !    cum = cum + ly(i + 1)
    !    if(p < 0) p = 0
    !    ly(i+1) = p       
    ! enddo

    ! cum = 0
    ! do i = 0, nz - 1
    !    p = min(ze, cum + lz(i+1) - 1) - max(zs, cum) + 1
    !    cum = cum + lz(i + 1)
    !    if(p < 0) p = 0
    !    lz(i+1) = p
    ! enddo

    call petsc_update_dist(lx, xs, xe)
    call petsc_update_dist(ly, ys, ye)
    call petsc_update_dist(lz, zs, ze)
    
    ! write(buffer, *) "lx=", lx, "ly=", ly, "lz=", lz, "\n"
    ! call PetscViewerASCIISynchronizedPrintf(&
    !      PETSC_VIEWER_STDOUT_WORLD, buffer, ierr)
    ! print*,"after:", "lx=", lx, "ly=",ly, "lz=", lz
    
    call DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
         DMDA_STENCIL_STAR, xe - xs + 1, ye - ys + 1, ze - zs + 1, &
         nx, ny, nz, 1, 0, lx, ly, lz, dst, ierr)

    CHKERRQ(ierr)

  end subroutine

  subroutine test_petsc_slice_dm()
    implicit none
#include "petsc.h"
    DM :: dst_dm, src_dm
    integer :: m, n, k
    integer :: ierr
    type(box_info) :: box
    
    m = 20
    n = 10
    k = 5
    
    call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,         &
         DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,                      &
         DMDA_STENCIL_STAR, m,n,k,PETSC_DECIDE,PETSC_DECIDE,     &
         PETSC_DECIDE, 1, 0,                                     &
         PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                  &
         PETSC_NULL_INTEGER, src_dm, ierr)

    box%rx = r(0, 11)
    box%ry = r(8, 9)
    box%rz = r(3, 4)
    
    call petsc_slice_dm(dst_dm, src_dm, box)

    box%rx = r(0, 5)
    box%ry = r(8, 9)
    box%rz = r(3, 4)
    
    call petsc_slice_dm(dst_dm, src_dm, box)
  end subroutine

  subroutine petsc_print(global_vec, prefix)
    implicit none
#include "petsc.h"
    Vec, intent(in) :: global_vec
    Vec :: local_vec
    VecScatter :: ctx
    integer :: ierr, v_shape(3)
    type(box_info) :: box
    PetscScalar, pointer :: arr(:)
    integer :: size, rank
    integer :: x, y, z, m, n, k
    character(len=*), optional :: prefix
    type(dist_info) :: dist
    DM :: global_dm
    PetscScalar,allocatable :: data3d(:,:,:)
    PetscScalar, pointer :: local_block_ptr(:)
    integer :: px, py, pz, accx, accy, accz, offset,block_size
    integer :: xs, xe, ys, ye, zs, ze, mm, nn, kk
    
    call VecScatterCreateToZero(global_vec, ctx, local_vec, ierr)
    call VecScatterBegin(ctx, global_vec, local_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr);
    call VecScatterEnd(ctx,global_vec,local_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr);

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    if(rank == 0) then
       call petsc_get_shape(v_shape, global_vec)
       call VecGetArrayF90(local_vec, arr, ierr)
       
       m = v_shape(1)
       n = v_shape(2)
       k = v_shape(3)
       
       allocate(data3d(m, n, k))
              
       call petsc_get_dm(global_dm, global_vec)
       call petsc_get_dist(dist, global_dm)

       offset = 0;
       accz = 0;
       do pz = 1, size(dist%lz)
          accy = 0
          do py = 1, size(dist%ly)
             accx = 0
             do px = 1, size(dist%lx)
                mm = dist%lx(px)
                nn = dist%ly(py)
                kk = dist%lz(pz)
                
                xs = accx+1; xe = xs + mm - 1 
                ys = accy+1; ye = ys + nn - 1
                zs = accz+1; ze = zs + kk - 1
                
                block_size = mm * nn * kk
                local_block_ptr &
                     => arr(offset+1:offset+block_size)

                ! print*, "m=", mm, "n=", nn, "k=",kk
                ! print*, "xs = ", xs, "ys=", ys, "zs=", zs
                ! print*, "xe = ", xe, "ye=", ye, "ze=", ze
                ! print*, "size(local_block_ptr) =", size(local_block_ptr)

                data3d(xs:xe,ys:ye,zs:ze) = &
                     reshape(local_block_ptr,(/mm, nn, kk/))
                
                offset = offset + block_size
                accx = accx + dist%lx(px)
             end do
             accy = accy + dist%ly(py)             
          end do
          accz = accz + dist%lz(pz)          
       end do
       
       if(present(prefix)) then
          write(*, "(4X, A)") prefix
       else
          write(*, "(4X, A)") "data = "
       endif

       do z = 1, k
          if(k > 1) &
               write(*,"(6X, A, I0.1)") "k = ", z
          do y = 1, n
             write(*,"(6X, 100g12.5)") &
                  data3d(:, y, z)
          end do
       end do
       
       ! do z = 1, k
       !    if(k > 1) &
       !         write(*,"(6X, A, I0.1)") "k = ", z
       !    do y = 1, n
       !       write(*,"(6X, 100g12.5)") &
       !            arr((z-1)*m*n+(y-1)*m+1:(z-1)*m*n+(y-1)*m+m)
       !    end do
       ! end do

       call VecRestoreArrayF90(local_vec, arr, ierr)
    end if

    call VecDestroy(local_vec, ierr)
    call VecScatterDestroy(ctx, ierr)
  end subroutine
end module

#undef I1
#undef I2
#undef I3
