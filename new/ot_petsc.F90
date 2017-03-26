
#define I1 :
#define I2 :,:
#define I3 :,:,:

module ot_petsc
  use ot_geom
  use ot_print
  
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
  end interface
  
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
  
!   subroutine dm_get_corners(DM_A, xs, ys, zs, xl, yl, zl)
!     implicit none
! #include "petsc.h"
!     DM :: DM_A
!     integer :: ierr
!     integer, intent(out) :: xs, ys, zs, xl, yl, zl
    
!     call DMDAGetCorners(DM_A, xs, ys, zs, xl, yl, zl,ierr)

!   end subroutine

  ! subroutine dm_get_corners(A, b)
  !   implicit none
  !   DM :: A
  !   integer :: xs, ys, zs, xl, yl, zl
  !   type(box_info), intent(out) :: b
  !   integer :: ierr
    
  !   call DMDAGetCorners(A, xs, ys, zs, xl, yl, zl,ierr)    
  !   b%starts = (/xs, ys, zs/)
  !   b%ends = (/xs+xl-1, ys+yl-1, zs+zl-1/)
  ! end subroutine

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

    b%starts = (/xs, ys, zs/)
    b%ends = (/xs+xl-1, ys+yl-1, zs+zl-1/)
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
    
    res%starts = 0
    res%ends = (/dx,dy,dz/) - 1
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
    call VecView(data, PETSC_VIEWER_STDOUT_WORLD,ierr)
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
    
    call DMGetGlobalVector(arr_dm, arr, ierr)
    call VecView(arr,PETSC_VIEWER_STDOUT_WORLD,ierr)    
  end subroutine

  subroutine petsc_new3d_by_dm(arr, arr_dm)
    implicit none
#include "petsc.h"    
    integer :: m, n, k, bm, bn, bk
    Vec, intent(out) :: arr
    DM, intent(in) :: arr_dm
    integer :: ierr
    
    call DMGetGlobalVector(arr_dm, arr, ierr)

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
    
    xs = box%starts(1)
    ys = box%starts(2)
    zs = box%starts(3)

    xe = box%ends(1)
    ye = box%ends(2)
    ze = box%ends(3)

    print*,"before:", "lx=", lx, "ly=",ly, "lz=", lz
    
    cum = 0
    do i = 0, nx - 1
       p = min(xe, cum + lx(i+1) - 1) - max(xs, cum) + 1
       cum = cum + lx(i + 1)
       if(p < 0) p = 0
       lx(i+1) = p
    enddo
    
    cum = 0
    do i = 0, ny - 1
       p = min(ye, cum + ly(i+1) - 1) - max(ys, cum) + 1
       cum = cum + ly(i + 1)
       if(p < 0) p = 0
       ly(i+1) = p       
    enddo

    cum = 0
    do i = 0, nz - 1
       p = min(ze, cum + lz(i+1) - 1) - max(zs, cum) + 1
       cum = cum + lz(i + 1)
       if(p < 0) p = 0
       lz(i+1) = p
    enddo

    ! write(buffer, *) "lx=", lx, "ly=", ly, "lz=", lz, "\n"
    ! call PetscViewerASCIISynchronizedPrintf(&
    !      PETSC_VIEWER_STDOUT_WORLD, buffer, ierr)
    print*,"after:", "lx=", lx, "ly=",ly, "lz=", lz
    
    call DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
         DMDA_STENCIL_STAR, xe - xs + 1, ye - ys + 1, ze - zs + 1, &
         nx, ny, nz, 1, 0, lx, ly, lz, dst, ierr)

    CHKERRQ(ierr)

    ! call PetscViewerASCIIPushSynchronized(PETSC_VIEWER_STDOUT_WORLD, ierr)

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

    box%starts = (/0, 8, 3/)
    box%ends   = (/11, 9, 4/)
    
    call petsc_slice_dm(dst_dm, src_dm, box)

    box%starts = (/0, 8, 3/)
    box%ends   = (/5, 9, 4/)
    
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
                  arr((z-1)*m*n+(y-1)*m+1:(z-1)*m*n+(y-1)*m+m)
          end do
       end do
       call VecRestoreArrayF90(local_vec, arr, ierr)
    end if

    call VecDestroy(local_vec, ierr)
    call VecScatterDestroy(ctx, ierr)
  end subroutine
end module

#undef I1
#undef I2
#undef I3
