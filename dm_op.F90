module dm_op
  use dm_type
  use dm
  implicit none
#include "mat_type.h"

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>

  
  type(Matrix)  :: MASK_X1, MASK_X2, MASK_Y1, MASK_Y2, MASK_Z1, MASK_Z2
  type(Matrix)  :: NAG_MASK_X1, NAG_MASK_X2, NAG_MASK_Y1, NAG_MASK_Y2, NAG_MASK_Z1, NAG_MASK_Z2
  type(Matrix)  :: REV_MASK_X1, REV_MASK_X2, REV_MASK_Y1, REV_MASK_Y2
  type(Matrix)  :: REV_MASK_Z1, REV_MASK_Z2
  type(Matrix)  :: HF_REV_MASK_Z1, HF_REV_MASK_Z2
  type(Matrix)  :: HF_MASK_Z1, HF_MASK_Z2
  type(Matrix)  :: ZEROS, ONES
  type(Matrix)  :: TRIL_IM_K1, TRIU_JM_K1
  
  private
  
  type(Matrix)  :: MAT_AXF, MAT_AXB, MAT_AYF, MAT_AYB, MAT_AZF, MAT_AZB
  type(Matrix)  :: MAT_DXF, MAT_DXB, MAT_DYF, MAT_DYB, MAT_DZF, MAT_DZB

  type(Matrix)  :: MAT_AXF_XY, MAT_AXB_XY, MAT_AYF_XY, MAT_AYB_XY
  type(Matrix)  :: MAT_DXF_XY, MAT_DXB_XY, MAT_DYF_XY, MAT_DYB_XY
  
  Mat  :: MAT_P, MAT_Q !used for shift a matrix in diagonal direction.
  Mat  :: MAT_R, MAT_T !used for shift a matrix in diagonal direction.
  Mat  :: UTI1, UTI2 !upper triangle I
  Mat  :: MAT_ALIGN_ROW
  Mat  :: MAT_SXD, MAT_SXU !left shift matrix, right shift matrix
  Mat  :: MAT_SYL, MAT_SYR !left shift matrix, right shift matrix
  
  public :: OP_AXF, OP_AXB, OP_AYF, OP_AYB
  public :: OP_DXF, OP_DXB, OP_DYF, OP_DYB
  
  public :: InitOperatorModule, FinalizeOperatorModule
  public :: AXF, AXB, AYF, AYB, AZF, AZB
  public :: DXF, DXB, DYF, DYB, DZF, DZB
  
  public :: DXC, DYC, DZC
  public :: CSUM , SHIFT

  public :: MASK_X1, MASK_X2, MASK_Y1, MASK_Y2, MASK_Z1, MASK_Z2
  public :: NAG_MASK_X1, NAG_MASK_X2, NAG_MASK_Y1, NAG_MASK_Y2, NAG_MASK_Z1, NAG_MASK_Z2
  public :: REV_MASK_X1, REV_MASK_X2, REV_MASK_Y1, REV_MASK_Y2
  public :: REV_MASK_Z1, REV_MASK_Z2
  public :: HF_REV_MASK_Z1, HF_REV_MASK_Z2
  public :: HF_MASK_Z1, HF_MASK_Z2
  public :: ZEROS, ONES
  public :: TRIL_IM_K1, TRIU_JM_K1
  
contains

  subroutine InitOperatorModule(m, n, k, isGlobal)
    implicit none

    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    logical :: is_global = .true.
    integer :: ierr, i
    integer, allocatable :: idxm(:), idxn(:), idxk(:)
    real(kind=8), allocatable :: zeros_mn(:), zeros_mk(:), zeros_nk(:)
    real(kind=8), allocatable :: ones_mn(:),  ones_mk(:),  ones_nk(:)
    integer :: COMM_TYPE
    integer :: ista, iend, num
    integer :: im, in, ik
    PetscScalar :: val
    
    allocate(idxm(m), idxn(n), idxk(k))
    allocate(ones_mn(m*n),  ones_mk(m*k),  ones_nk(n*k))
    allocate(zeros_mn(m*n), zeros_mk(m*k), zeros_nk(n*k))
    
    idxm = [(i, i=0,m-1)]
    idxn = [(i, i=0,n-1)]
    idxk = [(i, i=0,k-1)]

    ones_mn = 1
    ones_mk = 1
    ones_nk = 1

    zeros_mn = 0
    zeros_mk = 0
    zeros_nk = 0
    
    if(present(isGlobal)) is_global = isGlobal
    
    MAT_AXF = OP_AXF(m, n, k, is_global)
    MAT_AXB = OP_AXB(m, n, k, is_global)

    MAT_AXF_XY = OP_AXF(m, n, 1, is_global)
    MAT_AXB_XY = OP_AXB(m, n, 1, is_global)
    
    MAT_AYF = OP_AYF(m, n, k, is_global)
    MAT_AYB = OP_AYB(m, n, k, is_global)

    MAT_AYF_XY = OP_AYF(m, n, 1, is_global)
    MAT_AYB_XY = OP_AYB(m, n, 1, is_global)
    
    MAT_DXF = OP_DXF(m, n, k, is_global)
    MAT_DXB = OP_DXB(m, n, k, is_global)

    MAT_DXF_XY = OP_DXF(m, n, 1, is_global)
    MAT_DXB_XY = OP_DXB(m, n, 1, is_global)
    
    MAT_DYF = OP_DYF(m, n, k, is_global)
    MAT_DYB = OP_DYB(m, n, k, is_global)

    MAT_DYF_XY = OP_DYF(m, n, 1, is_global)
    MAT_DYB_XY = OP_DYB(m, n, 1, is_global)
    
    ONES = dm_ones(m, n, k, is_global)
    ZEROS = dm_zeros(m, n, k, is_global)
    
    MASK_X1 = ZEROS
    call dm_setvalues(MASK_X1, [(0)], idxn, idxk, ones_nk, ierr)

    MASK_X2 = ZEROS
    call dm_setvalues(MASK_X2, [(m-1)], idxn, idxk, ones_nk, ierr)

    REV_MASK_X1 = ONES
    call dm_setvalues(REV_MASK_X1, [(0)], idxn, idxk, zeros_nk, ierr)

    REV_MASK_X2 = ONES
    call dm_setvalues(REV_MASK_X2, [(m-1)], idxn, idxk, zeros_nk, ierr)
    
    MASK_Y1 = ZEROS
    call dm_setvalues(MASK_Y1, idxm, [(0)], idxk, ones_mk, ierr)

    MASK_Y2 = ZEROS
    call dm_setvalues(MASK_Y2, idxm, [(n-1)], idxk, ones_mk, ierr)
    
    REV_MASK_Y1 = ONES
    call dm_setvalues(REV_MASK_Y1, idxm, [(0)], idxk, zeros_mk, ierr)

    REV_MASK_Y2 = ONES
    call dm_setvalues(REV_MASK_Y2, idxm, [(n-1)], idxk, zeros_mk, ierr)
    
    MASK_Z1 = ZEROS
    call dm_setvalues(MASK_Z1, idxm, idxn, [(0)], ones_mn, ierr)

    MASK_Z2 = ZEROS
    call dm_setvalues(MASK_Z2, idxm, idxn, [(k-1)], ones_mn, ierr)
    
    REV_MASK_Z1 = ONES
    call dm_setvalues(REV_MASK_Z1, idxm, idxn, [(0)], zeros_mn, ierr)

    REV_MASK_Z2 = ONES
    call dm_setvalues(REV_MASK_Z2, idxm, idxn, [(k-1)], zeros_mn, ierr)

    NAG_MASK_X1 = - MASK_X1
    NAG_MASK_X2 = - MASK_X2
    NAG_MASK_Y1 = - MASK_Y1
    NAG_MASK_Y2 = - MASK_Y2
    NAG_MASK_Z1 = - MASK_Z1
    NAG_MASK_Z2 = - MASK_Z2
    
    HF_MASK_Z1 = 0.5 * MASK_Z1
    HF_MASK_Z2 = 0.5 * MASK_Z2

    HF_REV_MASK_Z1 = 0.5 * REV_MASK_Z1
    HF_REV_MASK_Z2 = 0.5 * REV_MASK_Z2
    
    if(is_global) then
       COMM_TYPE = PETSC_COMM_WORLD
    else
       COMM_TYPE = PETSC_COMM_SELF
    endif

    call mat_create(MAT_P, m, m, k, is_global, ierr)
    call mat_create(MAT_Q, n, n, k, is_global, ierr)
    
    val = 1
    call MatGetOwnershipRange(MAT_P, ista, iend, ierr)
    do im = ista, iend - 1
       ik = im / m
       if(im - m .ge. 0) then
         call MatSetValue(MAT_P, im, im - m, val, INSERT_VALUES, ierr)
       endif
    enddo

    call MatGetOwnershipRange(MAT_Q, ista, iend, ierr)
    do im = ista, iend - 1
       ik = im / n
       if(im + n .lt. n*k) then
          call MatSetValue(MAT_Q, im, im + n, val, INSERT_VALUES, ierr)
       endif
    enddo

    call mat_assemble(MAT_P, ierr)
    call mat_assemble(MAT_Q, ierr)

    call MatTranspose(MAT_P, MAT_INITIAL_MATRIX, MAT_R, ierr)
    call MatTranspose(MAT_Q, MAT_INITIAL_MATRIX, MAT_T, ierr)

    call mat_assemble(MAT_R, ierr)
    call mat_assemble(MAT_T, ierr)

    ! call mat_view(MAT_P, ierr)
    ! call mat_view(MAT_Q, ierr)
    ! call mat_view(MAT_R, ierr)
    ! call mat_view(MAT_T, ierr)
    
    !**********************************************************
    ! [I I I 0]   [A 0 0 0]   [A B C 0]
    ! [0 I I 0] * [0 B 0 0] = [0 B C 0]
    ! [0 0 I 0]   [0 0 C 0]   [0 0 C 0]
    ! [0 0 0 0]   [0 0 0 D]   [0 0 0 0]
    !**********************************************************    
    call mat_create(UTI1, m, m, k, is_global, ierr)
    call MatGetOwnershipRange(UTI1, ista, iend, ierr)
    do im = ista, iend-1
       ik = im / m
       do in = 0,k-ik-2
          call MatSetValue(UTI1, im, im + in * m, val, INSERT_VALUES, ierr)
       enddo
    enddo
    call mat_assemble(UTI1, ierr)
    ! call mat_view(UTI1, ierr)
    
    !**********************************************************
    ! [0 I I I]   [A 0 0 0]   [0 B C D]
    ! [0 0 I I] * [0 B 0 0] = [0 0 C D]
    ! [0 0 0 I]   [0 0 C 0]   [0 0 0 D]
    ! [0 0 0 0]   [0 0 0 D]   [0 0 0 0]
    !**********************************************************    
    call mat_create(UTI2, m, m, k, is_global, ierr)
    call MatGetOwnershipRange(UTI2, ista, iend, ierr)
    do im = ista, iend-1
       ik = im / m
       do in = 1,k-ik-1
          call MatSetValue(UTI2, im, im + in * m, val, INSERT_VALUES, ierr)
       enddo
    enddo
    call mat_assemble(UTI2, ierr)
    ! call mat_view(UTI2, ierr)
    
    !**********************************************************    
    ! [I I I I]   [A 0 0 0]   [A B C D]
    ! [0 0 0 0] * [0 B 0 0] = [0 0 0 0]
    ! [0 0 0 0]   [0 0 C 0]   [0 0 0 0]
    ! [0 0 0 0]   [0 0 0 D]   [0 0 0 0]
    !**********************************************************    
    call mat_create(MAT_ALIGN_ROW, m, m, k, is_global, ierr)
    call MatGetOwnershipRange(MAT_ALIGN_ROW, ista, iend, ierr)
    do im = ista, iend-1
       ik = im / m
       if(ik .gt. 0) exit
       do in = 0,k-ik-1
          call MatSetValue(MAT_ALIGN_ROW, im, im + in*m, &
               val, INSERT_VALUES, ierr)
       enddo
    enddo
    call mat_assemble(MAT_ALIGN_ROW, ierr)
    !call mat_view(MAT_ALIGN_ROW, ierr)

    !***********************************************************
    ! [0 1 0 0]
    ! [0 0 1 0] 
    ! [0 0 0 1]
    ! [0 0 0 0]
    ! MAT_SXD * A ==> Shift downward
    !***********************************************************
    call mat_create(MAT_SXD, m, m, k, is_global, ierr)
    val = 1.0
    call MatGetOwnershipRange(MAT_SXD, ista, iend, ierr)
    do im = ista, iend - 1
       ik = im / m
       if(mod(im, m) .lt. m-1) then
          call MatSetValue(MAT_SXD, im, mod(im,m) + 1 + ik * m, val, &
               INSERT_VALUES, ierr)
       endif
    enddo
    call mat_assemble(MAT_SXD, ierr)
    !call mat_view(MAT_SXD, ierr)
    
    !***********************************************************
    ! [0 0 0 0]
    ! [1 0 0 0] 
    ! [0 1 0 0]
    ! [0 0 1 0]
    ! MAT_SXU * A ==> Shift upward
    !***********************************************************
    call mat_create(MAT_SXU, m, m, k, is_global, ierr)
    val = 1.0
    call MatGetOwnershipRange(MAT_SXU, ista, iend, ierr)
    do im = ista, iend - 1
       ik = im / m
       if(mod(im, m) .gt. 0) then
          call MatSetValue(MAT_SXU, im, mod(im,m) - 1 + ik * m, val, &
               INSERT_VALUES, ierr)
       endif
    enddo
    call mat_assemble(MAT_SXU, ierr)
    !call mat_view(MAT_SXU, ierr)

    !***********************************************************
    ! [0 1 0 0]
    ! [0 0 1 0] 
    ! [0 0 0 1]
    ! [0 0 0 0]
    ! A * MAT_SYR ==> Shift to right
    !***********************************************************
    call mat_create(MAT_SYR, n, n, k, is_global, ierr)
    val = 1.0
    call MatGetOwnershipRange(MAT_SYR, ista, iend, ierr)
    do im = ista, iend - 1
       ik = im / n
       if(mod(im, n) .lt. n-1) then
          call MatSetValue(MAT_SYR, im, mod(im,n) + 1 + ik * n, val, &
               INSERT_VALUES, ierr)
       endif
    enddo
    call mat_assemble(MAT_SYR, ierr)
    !call mat_view(MAT_SYR, ierr)

    !***********************************************************
    ! [0 0 0 0]
    ! [1 0 0 0] 
    ! [0 1 0 0]
    ! [0 0 1 0]
    ! A * MAT_SYL ==> Shift to left
    !***********************************************************
    call mat_create(MAT_SYL, n, n, k, is_global, ierr)
    val = 1.0
    call MatGetOwnershipRange(MAT_SYL, ista, iend, ierr)
    do im = ista, iend - 1
       ik = im / n
       if(mod(im, n) .gt. 0) then
          call MatSetValue(MAT_SYL, im, mod(im,n) - 1 + ik * n, val, &
               INSERT_VALUES, ierr)
       endif
    enddo
    call mat_assemble(MAT_SYL, ierr)
    !call mat_view(MAT_SYL, ierr)

    !***********************************************************
    ! [1 0 0 0]
    ! [1 1 0 0] 
    ! [1 1 1 0]
    ! [1 1 1 1]
    !***********************************************************    
    call dm_create(TRIL_IM_K1, m, m, 1, is_global, ierr)
    val = 1.0
    call MatGetOwnershipRange(TRIL_IM_K1%x, ista, iend, ierr)
    do im = ista, iend - 1
       ik = im / m
       do in = 0, mod(im, m)
          if(in .le. mod(im, m)) then
             call MatSetValue(TRIL_IM_K1%x, im, in + ik * m, val, &
                  INSERT_VALUES, ierr)
          endif
       enddo
    enddo
    call mat_assemble(TRIL_IM_K1%x, ierr)
    !call mat_view(TRIL_IM_K1%x, ierr)

    !***********************************************************
    ! [1 1 1 1]
    ! [0 1 1 1] 
    ! [0 1 1 1]
    ! [0 0 0 1]
    !***********************************************************    
    call dm_create(TRIU_JM_K1, n, n, 1, is_global, ierr)
    val = 1.0
    call MatGetOwnershipRange(TRIU_JM_K1%x, ista, iend, ierr)
    do im = ista, iend - 1
       ik = im / n
       do in = 0, mod(im, n)
          if(in .le. mod(im, n)) then
             call MatSetValue(TRIU_JM_K1%x, im, in + ik * n, val, &
                  INSERT_VALUES, ierr)
          endif
       enddo
    enddo
    call mat_assemble(TRIU_JM_K1%x, ierr)
    !call dm_view(TRIU_JM_K1, ierr)
    
  end subroutine 

  subroutine FinalizeOperatorModule()
    integer :: ierr

    call dm_destroy(MAT_AXF, ierr)
    call dm_destroy(MAT_AXB, ierr)
    call dm_destroy(MAT_AYF, ierr)
    call dm_destroy(MAT_AYB, ierr)
    call dm_destroy(MAT_DXF, ierr)
    call dm_destroy(MAT_DXB, ierr)
    call dm_destroy(MAT_DYF, ierr)
    call dm_destroy(MAT_DYB, ierr)
    call dm_destroy(MASK_X1, ierr)
    call dm_destroy(MASK_X2, ierr)
    call dm_destroy(MASK_Y1, ierr)
    call dm_destroy(MASK_Y2, ierr)
    call dm_destroy(MASK_Z1, ierr)
    call dm_destroy(MASK_Z2, ierr)

    call dm_destroy(MAT_AXF_XY, ierr)
    call dm_destroy(MAT_AXB_XY, ierr)
    call dm_destroy(MAT_AYF_XY, ierr)
    call dm_destroy(MAT_AYB_XY, ierr)
    call dm_destroy(MAT_DXF_XY, ierr)
    call dm_destroy(MAT_DXB_XY, ierr)
    call dm_destroy(MAT_DYF_XY, ierr)
    call dm_destroy(MAT_DYB_XY, ierr)

    call dm_destroy(NAG_MASK_X1, ierr)
    call dm_destroy(NAG_MASK_X2, ierr)
    call dm_destroy(NAG_MASK_Y1, ierr)
    call dm_destroy(NAG_MASK_Y2, ierr)
    call dm_destroy(NAG_MASK_Z1, ierr)
    call dm_destroy(NAG_MASK_Z2, ierr)
        
    call dm_destroy(REV_MASK_X1, ierr)
    call dm_destroy(REV_MASK_X2, ierr)
    call dm_destroy(REV_MASK_Y1, ierr)
    call dm_destroy(REV_MASK_Y2, ierr)
    call dm_destroy(REV_MASK_Z1, ierr)
    call dm_destroy(REV_MASK_Z2, ierr)

    call dm_destroy(HF_MASK_Z1, ierr)
    call dm_destroy(HF_MASK_Z2, ierr)
    call dm_destroy(HF_REV_MASK_Z1, ierr)
    call dm_destroy(HF_REV_MASK_Z2, ierr)
    
    call dm_destroy(ZEROS, ierr)
    call dm_destroy(ONES, ierr)

    call mat_destroy(MAT_P, ierr)
    call mat_destroy(MAT_Q, ierr)
    call mat_destroy(MAT_R, ierr)
    call mat_destroy(MAT_T, ierr)

    call mat_destroy(UTI1, ierr)
    call mat_destroy(UTI2, ierr)

    call mat_destroy(MAT_SXD, ierr)
    call mat_destroy(MAT_SXU, ierr)
    call mat_destroy(MAT_SYL, ierr)
    call mat_destroy(MAT_SYR, ierr)

    call dm_destroy(TRIL_IM_K1, ierr)
    call dm_destroy(TRIU_JM_K1, ierr)
    
    call mat_destroy(MAT_ALIGN_ROW, ierr)

  end subroutine 

  !****************************************
  !          Oparator matrices
  !***************************************
  
  ! OP_AXF= 0.5*[1  1  0  0  0  0  0]   
  !             [0  1  1  0  0  0  0]   
  !             [0  0  1  1  0  0  0]   
  !             [0  0  0  1  1  0  0]   
  !             [0  0  0  0  1  1  0]   
  !             [0  0  0  0  0  1  1]   
  !             [0  0  0  0  0  0  2]
  function OP_AXF(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    
    res = dm_eye(m, m, k, is_global)
    L2 = dm_zeros(m, 1, k, is_global) .yj. dm_eye(m, m-1, k, is_global)

    call dm_setvalues(res, [(m-1)], [(m-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    res = (res + L2) * .5
    
    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
  end function 

  ! OP_AYF= 0.5*[1  0  0  0  0]        
  !             [1  1  0  0  0]        
  !             [0  1  1  0  0]        
  !             [0  0  1  1  0]        
  !             [0  0  0  1  2]        
  ! Y*OP_AYF= 0.5*[ Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15  2Y15 ]
  !               [ Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25  2Y25 ]
  !               [ Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35  2Y35 ]
  !               [ Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45  2Y45 ]
  !               [ Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55  2Y55 ]
  !               [ Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65  2Y65 ]
  !               [ Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75  2Y75 ]
  function OP_AYF(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    res = dm_eye(n, n,   k, is_global)
    L2 = dm_zeros(1, n, k, is_global) .xj. dm_eye(n-1, n, k, is_global)

    call dm_setvalues(res, [(n-1)], [(n-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    res = (res + L2) * 0.5

    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
  end function 

  !  OP_AXB= 0.5*[2  0  0  0  0  0  0]   
  !              [1  1  0  0  0  0  0]   
  !              [0  1  1  0  0  0  0]   
  !              [0  0  1  1  0  0  0]   
  !              [0  0  0  1  1  0  0]   
  !              [0  0  0  0  1  1  0]   
  !              [0  0  0  0  0  1  1]   
  function OP_AXB(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    
    res = dm_eye(m,  m,   k, is_global)
    L2 = dm_zeros(1, m, k, is_global) .xj. dm_eye(m-1, m, k, is_global)

    call dm_setvalues(res, [(0)], [(0)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    res = (res + L2) * 0.5

    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
  end function OP_AXB

  ! OP_AYB= 0.5*[2  1  0  0  0]        
  !             [0  1  1  0  0]        
  !             [0  0  1  1  0]        
  !             [0  0  0  1  1]        
  !             [0  0  0  0  1]        
  !                                    
  ! Y*OP_AYB= 0.5*[ 2Y11  Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15 ]
  !             [ 2Y21  Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25 ]
  !             [ 2Y31  Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35 ]
  !             [ 2Y41  Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45 ]
  !             [ 2Y51  Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55 ]
  !             [ 2Y61  Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65 ]
  !             [ 2Y71  Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75 ]
  function OP_AYB(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal

    res = dm_eye(n, n, k, is_global)
    L2  = dm_zeros(n, 1, k, is_global) .yj. dm_eye(n, n-1, k, is_global)

    call dm_setvalues(res, [(0)], [(0)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    res = (res + L2) * 0.5
    
    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)

  end function OP_AYB


  !OP_DXB1_XY=[  0  0  0  0  0  0  0]     
  !           [ -1  1  0  0  0  0  0]      
  !           [  0 -1  1  0  0  0  0]      
  !           [  0  0 -1  1  0  0  0]      
  !           [  0  0  0 -1  1  0  0]      
  !           [  0  0  0  0 -1  1  0]      
  !           [  0  0  0  0  0 -1  1]
  ! OP_DXB1_XY*X =[       0        0         0         0         0]
  !               [ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
  !               [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
  !               [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
  !               [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]  
  !               [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
  !               [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]
  function OP_DXB(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    res = dm_eye(m, m, k, is_global)
    L2  = dm_zeros(1, m, k, is_global) .xj. (-1.0*dm_eye(m-1, m, k, is_global))

    call dm_setvalues(res, [(0)], [(0)], [(i, i=0,k-1)], [(0,i=0,k-1)], ierr)

    res = res + L2
    ! call dm_destroy(L1, ierr)
    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
    
  end function OP_DXB

  ! OP_DYB1_XY=[  0 -1  0  0  0 ]         
  !            [  0  1 -1  0  0 ]         
  !            [  0  0  1 -1  0 ]         
  !            [  0  0  0  1 -1 ]         
  !            [  0  0  0  0  1 ]         
  ! Y*OP_DYB1_XY =[ 0  Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14 ]
  !               [ 0  Y22-Y21   Y23-Y22   Y24-Y23   Y15-Y14 ]
  !               [ 0  Y32-Y31   Y33-Y32   Y34-Y33   Y35-Y34 ]
  !               [ 0  Y42-Y41   Y43-Y42   Y44-Y43   Y45-Y44 ]
  !               [ 0  Y52-Y51   Y53-Y52   Y54-Y53   Y55-Y54 ]
  !               [ 0  Y62-Y61   Y63-Y62   Y64-Y63   Y65-Y64 ]
  !               [ 0  Y72-Y71   Y73-Y72   Y74-Y73   Y75-Y74 ]
  function OP_DYB(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    
    res = dm_eye(n, n, k, is_global)
    L2 = dm_zeros(n, 1, k, is_global) .yj. (-1.0*dm_eye(n, n-1, k, is_global))

    call dm_setvalues(res, [(0)], [(0)], [(i, i=0,k-1)], [(0,i=0,k-1)], ierr)

    res = res + L2

    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
    
  end function OP_DYB

  ! OP_DXF1_XY=[ -1  1  0  0  0  0  0]    
  !            [  0 -1  1  0  0  0  0]    
  !            [  0  0 -1  1  0  0  0]    
  !            [  0  0  0 -1  1  0  0]    
  !            [  0  0  0  0 -1  1  0]    
  !            [  0  0  0  0  0 -1  1]    
  !            [  0  0  0  0  0  0  0]
  ! OP_DXF1_XY*X =[ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
  !               [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
  !               [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
  !               [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
  !               [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
  !               [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]  
  !               [       0        0         0         0         0 ]
  function OP_DXF(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    res = -1.0 * dm_eye(m, m, k, is_global)
    L2 = dm_zeros(m, 1, k, is_global) .yj. dm_eye(m, m-1, k, is_global)

    res = res + L2

    call dm_setvalues(res, [(m-1)], [(m-1)], [(i, i=0,k-1)], [(0,i=0,k-1)], ierr)
    
    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
    
  end function OP_DXF
  
  ! OP_DYF1_XY=[ -1  0  0  0  0 ]        
  !            [  1 -1  0  0  0 ]        
  !            [  0  1 -1  0  0 ]        ]
  !            [  0  0  1 -1  0 ]        
  !            [  0  0  0  1  0 ]        
  !                                      
  !                                      
  ! Y*OP_DYF1_XY =[ Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14  0 ]
  !               [ Y22-Y21   Y23-Y22   Y24-Y23   Y25-Y24  0 ]
  !               [ Y22-Y21   Y33-Y32   Y34-Y33   Y35-Y34  0 ]
  !               [ Y22-Y21   Y43-Y42   Y44-Y43   Y45-Y44  0 ]
  !               [ Y22-Y21   Y53-Y52   Y54-Y53   Y55-Y54  0 ]
  !               [ Y22-Y21   Y63-Y62   Y64-Y63   Y65-Y64  0 ]  
  !               [ Y22-Y21   Y73-Y72   Y74-Y73   Y75-Y74  0 ]
  function OP_DYF(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    res = -1.0 * dm_eye(n, n, k, is_global)
    L2 = dm_zeros(1, n, k, is_global) .xj. dm_eye(n-1, n, k, is_global)

    call dm_setvalues(res, [(n-1)], [(n-1)], [(i, i=0,k-1)], [(0,i=0,k-1)], ierr)

    res = res + L2
    ! call dm_destroy(L1, ierr)
    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
  end function OP_DYF

  !****************************************
  !            Oparations
  !***************************************
  function AXF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    if(A%nz /= 1) then
       res = MAT_AXF * A
    else
       res = MAT_AXF_XY * A       
    endif

    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function AXB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    if(A%nz /= 1) then
       res = MAT_AXB * A
    else
       res = MAT_AXB_XY * A       
    endif
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function AYF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    if(A%nz /= 1) then
       res = A * MAT_AYF
    else
       res = A * MAT_AYF_XY
    endif
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 
  
  function AYB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    if(A%nz /= 1) then
       res = A * MAT_AYB
    else
       res = A * MAT_AYB_XY
    endif
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function AZF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    call mat_assemble(A%x, ierr)
    
    call MatMatMatMult(MAT_R, A%x,  MAT_T, MAT_INITIAL_MATRIX, &
         PETSC_DEFAULT_REAL, res%x, ierr)
    
    res = 0.5*((A+res)+(A .em. MASK_Z2))

    res%nx = A%nx
    res%ny = A%ny
    res%nz = A%nz
    res%isGlobal = A%isGlobal
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function AZB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    call mat_assemble(A%x, ierr)
    
    call MatMatMatMult(MAT_P, A%x,  MAT_Q, MAT_INITIAL_MATRIX, &
         PETSC_DEFAULT_REAL, res%x, ierr)
    
    res = 0.5*((A+res) +(A .em. MASK_Z1))

    res%nx = A%nx
    res%ny = A%ny
    res%nz = A%nz
    res%isGlobal = A%isGlobal
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)    
  end function 
  
  function DXF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    if(A%nz /= 1) then
       res = MAT_DXF * A
    else
       res = MAT_DXF_XY * A
    endif
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function DXB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    if(A%nz /= 1) then
       res = MAT_DXB * A
    else
       res = MAT_DXB_XY * A
    endif
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function DXC(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = MAT_DXB * MAT_AXF * A
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 
  
  function DYF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    if(A%nz /= 1) then
       res = A * MAT_DYF
    else
       res = A * MAT_DYF_XY
    endif

    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function DYB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    if(A%nz /= 1) then
       res = A * MAT_DYB
    else
       res = A * MAT_DYB_XY
    endif
    
    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function DYC(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = A * MAT_AYF * MAT_DYB

    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function DZF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    call mat_assemble(A%x, ierr)
    
    call MatMatMatMult(MAT_R, A%x,  MAT_T, MAT_INITIAL_MATRIX, &
         PETSC_DEFAULT_REAL, res%x, ierr)

    res%nx = A%nx
    res%ny = A%ny
    res%nz = A%nz
    res%isGlobal = A%isGlobal
    
    res = (res - A) .em. REV_MASK_Z2

    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function 

  function DZB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr

    call mat_assemble(A%x, ierr)
    
    call MatMatMatMult(MAT_P, A%x,  MAT_Q, MAT_INITIAL_MATRIX, &
         PETSC_DEFAULT_REAL, res%x, ierr)

    res%nx = A%nx
    res%ny = A%ny
    res%nz = A%nz
    res%isGlobal = A%isGlobal
    
    res = (A-res) .em. REV_MASK_Z1

    ! if(A%xtype == MAT_XTYPE_IMPLICIT) then
    !    call dm_destroy(A, ierr)
    ! endif
    
    call dm_set_implicit(res, ierr)
  end function 

  function DZC(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res, W
    integer :: ierr

    call mat_assemble(A%x, ierr)

    call MatMatMatMult(MAT_P, A%x, MAT_Q, MAT_INITIAL_MATRIX, &
         PETSC_DEFAULT_REAL, res%x, ierr)

    call MatMatMatMult(MAT_R, A%x, MAT_T, MAT_INITIAL_MATRIX, &
         PETSC_DEFAULT_REAL, W%x, ierr)

    W%nx = A%nx
    W%ny = A%ny
    W%nz = A%nz
    W%isGlobal = A%isGlobal
    
    res%nx = A%nx
    res%ny = A%ny
    res%nz = A%nz
    res%isGlobal = A%isGlobal
    
    res = W - res

    res = (res .em. HF_REV_MASK_Z1) + (A .em. HF_MASK_Z2)

    ! if(A%xtype == MAT_XTYPE_IMPLICIT) then
    !    call dm_destroy(A, ierr)
    ! endif
    
    call dm_set_implicit(res, ierr)
  end function 

  function SHIFT(A, axis, direction) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer, intent(in) :: axis
    integer, intent(in) :: direction
    integer :: ierr

    res%nx = A%nx
    res%ny = A%ny
    res%nz = A%nz
    res%isGlobal = A%isGlobal
    
    call mat_assemble(A%x, ierr)
    
    if(axis .eq. 1) then
       if(direction .eq. 1) then
          call MatMatMult(MAT_SXD, A%x, MAT_INITIAL_MATRIX, &
               PETSC_DEFAULT_REAL, res%x, ierr)
       else if(direction .eq. -1) then
          call MatMatMult(MAT_SXU, A%x, MAT_INITIAL_MATRIX, &
               PETSC_DEFAULT_REAL, res%x, ierr)
       endif
    else if(axis .eq. 2) then
       if(direction .eq. 1) then
          call MatMatMult(A%x, MAT_SYR, MAT_INITIAL_MATRIX, &
               PETSC_DEFAULT_REAL, res%x, ierr)
       else if(direction .eq. -1) then
          call MatMatMult(A%x, MAT_SYL, MAT_INITIAL_MATRIX, &
               PETSC_DEFAULT_REAL, res%x, ierr)
       endif
    else
       if(direction .eq. 1) then
          call MatMatMatMult(MAT_R, A%x, MAT_T, MAT_INITIAL_MATRIX, &
               PETSC_DEFAULT_REAL, res%x, ierr)
       else if(direction .eq. -1) then
          call MatMatMatMult(MAT_P, A%x, MAT_Q, MAT_INITIAL_MATRIX, &
               PETSC_DEFAULT_REAL, res%x, ierr)
       endif
    endif

    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif

    call dm_set_implicit(res, ierr)
  end function
  
  !< Cumulative summation
  function CSUM(A, type) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: type
    Mat :: W
    integer      :: ierr
    integer :: im, in, ik, ista, iend, j, in1
    integer :: from, to, count
    integer :: m, n, k, col
    PetscInt, allocatable :: idxn(:)
    PetscScalar, allocatable :: row(:)
    logical :: isGlobal
    !res = ZEROS
    integer :: COMM_TYPE

    m = A%nx
    n = A%ny
    k = A%nz

    call mat_assemble(A%x, ierr)

    if(type == 1) then
       call MatMatMult(UTI1, A%x, MAT_INITIAL_MATRIX, &
            PETSC_DEFAULT_REAL, W, ierr)
    else if(type == 2) then
       call MatMatMult(UTI2, A%x, MAT_INITIAL_MATRIX, &
            PETSC_DEFAULT_REAL, W, ierr)
    else if(type == 3 .or. type == 4) then
       call MatMatMult(MAT_ALIGN_ROW, A%x, MAT_INITIAL_MATRIX, &
            PETSC_DEFAULT_REAL, W, ierr)
    else
       call dm_printf("Error: Unknown sum type.", ierr)
       stop
    endif
    
    call MatGetOwnershipRange(W, ista, iend, ierr)

    allocate(row(A%ny * A%nz), idxn(n*k))

    if(type == 4) then
       call dm_create(res, A%nx, A%ny, 1, A%isGlobal, ierr) 
    else
       call dm_create(res, A%nx, A%ny, A%nz, A%isGlobal, ierr)
    endif

    idxn = 0
    
    do im = ista, iend-1
       ik = im / m

       if(type == 3 .and. ik > 0) exit
       if(type == 4 .and. ik > 0) exit
       
       call MatGetRow(W, im, col, idxn, row, ierr)

       from = 1
       count = 1
       
       do j = 2,col
          in= idxn(j)
          in1 = idxn(j-1)
          if( (in/n .ne. in1/n) .or. (j == col)) then

             if(j==col) count = count + 1
             to = from + count - 1
             
             call MatSetValues(res%x, 1, im, count, &
                  ik*n + mod(idxn(from:to), n), &
                  row(from:to), ADD_VALUES, ierr)

             from = from + count
             count = 0
          endif
          count = count + 1
       enddo
       
       call MatRestoreRow(W, im, col, idxn, row, ierr)
    enddo
    
    deallocate(row, idxn)
    call mat_destroy(W, ierr)

    if(A%xtype == MAT_XTYPE_IMPLICIT) then
       call dm_destroy(A, ierr)
    endif
    call dm_set_implicit(res, ierr)
  end function 
  
  end module dm_op
