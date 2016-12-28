module dm_op
  use dm_type
  use dm
  implicit none
#include "mat_type.h"

  private
  
  type(Matrix)  :: MAT_AXF, MAT_AXB, MAT_AYF, MAT_AYB, MAT_AZF, MAT_AZB
  type(Matrix)  :: MAT_DXF, MAT_DXB, MAT_DYF, MAT_DYB, MAT_DZF, MAT_DZB

  type(Matrix)  :: MASK_X1, MASK_X2, MASK_Y1, MASK_Y2, MASK_Z1, MASK_Z2
  type(Matrix)  :: REV_MASK_X1, REV_MASK_X2, REV_MASK_Y1, REV_MASK_Y2
  type(Matrix)  :: REV_MASK_Z1, REV_MASK_Z2
  type(Matrix)  :: ZEROS, ONES
  
  public :: OP_AXF, OP_AXB, OP_AYF, OP_AYB, OP_AZF, OP_AZB
  public :: OP_DXF, OP_DXB, OP_DYF, OP_DYB, OP_DZF, OP_DZB
  
  public :: CreateOperators
  public :: AXF, AXB, AYF, AYB, AZF, AZB
  public :: DXF, DXB, DYF, DYB, DZF, DZB
  public :: DXC, DYC, DZC

contains

  subroutine CreateOperators(m, n, k, isGlobal)

    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    logical :: is_global = .true.
    integer :: ierr, i
    
    if(present(isGlobal)) is_global = isGlobal
    
    MAT_AXF = OP_AXF(m, n, k, is_global)
    MAT_AXB = OP_AXB(m, n, k, is_global)
    
    MAT_AYF = OP_AYF(m, n, k, is_global)
    MAT_AYB = OP_AYB(m, n, k, is_global)
    
    MAT_AZF = OP_AZF(m, n, k, is_global)
    MAT_AZB = OP_AZB(m, n, k, is_global)
    
    MAT_DXF = OP_DXF(m, n, k, is_global)
    MAT_DXB = OP_DXB(m, n, k, is_global)

    MAT_DYF = OP_DYF(m, n, k, is_global)
    MAT_DYB = OP_DYB(m, n, k, is_global)

    MAT_DZF = OP_DZF(m, n, k, is_global)
    MAT_DZB = OP_DZB(m, n, k, is_global)

    ONES = dm_ones(m, n, k, is_global)
    ZEROS = dm_zeros(m, n, k, is_global)
    
    MASK_X1 = ZEROS
    call dm_setvalues(MASK_X1, [(0)], [(i,i=0,n-1)], [(i,i=0,k-1)], [(1,i=0,k*n-1)], ierr)

    MASK_X2 = ZEROS
    call dm_setvalues(MASK_X2, [(m-1)], [(i,i=0,n-1)], [(i,i=0,k-1)], [(1,i=0,k*n-1)], ierr)

    REV_MASK_X1 = ONES
    call dm_setvalues(REV_MASK_X1, [(0)], [(i,i=0,n-1)], [(i,i=0,k-1)], [(0,i=0,k*n-1)], ierr)

    REV_MASK_X2 = ONES
    call dm_setvalues(REV_MASK_X2, [(m-1)], [(i,i=0,n-1)], [(i,i=0,k-1)], [(0,i=0,k*n-1)], ierr)
    
    MASK_Y1 = ZEROS
    call dm_setvalues(MASK_Y1, [(i,i=0,m-1)], [(0)], [(i,i=0,k-1)], [(1,i=0,m*k-1)], ierr)

    MASK_Y2 = ZEROS
    call dm_setvalues(MASK_Y2, [(i,i=0,m-1)], [(n-1)], [(i,i=0,k-1)], [(1,i=0,m*k-1)], ierr)
    
    REV_MASK_Y1 = ONES
    call dm_setvalues(REV_MASK_Y1, [(i,i=0,m-1)], [(0)], [(i,i=0,k-1)], [(0,i=0,m*k-1)], ierr)

    REV_MASK_Y2 = ONES
    call dm_setvalues(REV_MASK_Y2, [(i,i=0,m-1)], [(n-1)], [(i,i=0,k-1)], [(0,i=0,m*k-1)], ierr)
    
    MASK_Z1 = ZEROS
    call dm_setvalues(MASK_Z1, [(i,i=0,m-1)], [(i,i=0,n-1)], [(0)], [(1,i=0,m*n-1)], ierr)

    MASK_Z2 = ZEROS
    call dm_setvalues(MASK_Z2, [(i,i=0,m-1)], [(i,i=0,n-1)], [(k-1)], [(1,i=0,m*n-1)], ierr)
    
    REV_MASK_Z1 = ONES
    call dm_setvalues(REV_MASK_Z1, [(i,i=0,m-1)], [(i,i=0,n-1)], [(0)], [(0,i=0,m*n-1)], ierr)

    REV_MASK_Z2 = ONES
    call dm_setvalues(REV_MASK_Z2, [(i,i=0,m-1)], [(i,i=0,n-1)], [(k-1)], [(0,i=0,m*n-1)], ierr)
    
  end subroutine 
  
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
  !                                    
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
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    res = dm_eye(n, n,   k, is_global)
    L2 = dm_zeros(1, n, k, is_global) .xj. dm_eye(n-1, n, k, is_global)

    call dm_setvalues(res, [(n-1)], [(n-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    res = (res + L2) * 0.5

    !call dm_destroy(L1, ierr)
    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
  end function 

  function OP_AZF_P(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    ! L1 = dm_eye(m, m, k, is_global)
    ! L2 = dm_zeros(m, 1, k, is_global) .yj. dm_eye(m, n-1, k, is_global)

    ! call dm_setvalues(L1, [(m-1)], [(n-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    ! res = L1 + L2
    ! call dm_destroy(L1, ierr)
    ! call dm_destroy(L2, ierr)

    call dm_set_implicit(res, ierr)
  end function OP_AZF_P

  function OP_AZF(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    ! L1 = dm_eye(m, n, k, is_global)
    ! L2 = dm_zeros(m, 1, k, is_global) .yj. dm_eye(m, n-1, k, is_global)

    ! call dm_setvalues(L1, [(m-1)], [(n-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    ! res = L1 + L2
    ! call dm_destroy(L1, ierr)
    ! call dm_destroy(L2, ierr)

    call dm_set_implicit(res, ierr)
  end function OP_AZF
  

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
    type(Matrix) :: L1, L2
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
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal

    res = dm_eye(n, n, k, is_global)
    L2  = dm_zeros(n, 1, k, is_global) .yj. dm_eye(n, n-1, k, is_global)

    call dm_setvalues(res, [(0)], [(0)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    res = (res + L2) * 0.5
    !call dm_destroy(L1, ierr)
    
    call dm_destroy(L2, ierr)
    call dm_set_implicit(res, ierr)
  end function OP_AYB


  function OP_AZB(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    ! L1 = dm_eye(m, n, k, is_global)
    ! L2 = dm_zeros(m, 1, k, is_global) .yj. dm_eye(m, n-1, k, is_global)

    ! call dm_setvalues(L1, [(m-1)], [(n-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    ! res = L1 + L2
    ! call dm_destroy(L1, ierr)
    ! call dm_destroy(L2, ierr)

    call dm_set_implicit(res, ierr)
  end function OP_AZB


  !OP_DXB1_XY=[  0  0  0  0  0  0  0]     
  !           [ -1  1  0  0  0  0  0]      
  !           [  0 -1  1  0  0  0  0]      
  !           [  0  0 -1  1  0  0  0]      
  !           [  0  0  0 -1  1  0  0]      
  !           [  0  0  0  0 -1  1  0]      
  !           [  0  0  0  0  0 -1  1]

  ! OP_DXB1_XY*X =[       0        0         0         0         0]
  !             [ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
  !             [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
  !             [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
  !             [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]  
  !             [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
  !             [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]
  
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


! % OP_DYB1_XY=[  0 -1  0  0  0 ]         
! %            [  0  1 -1  0  0 ]         
! %            [  0  0  1 -1  0 ]         
! %            [  0  0  0  1 -1 ]         
! %            [  0  0  0  0  1 ]         
! %                                       

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


  function OP_DZB(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    ! L1 = dm_eye(m, n, k, is_global)
    ! L2 = dm_zeros(m, 1, k, is_global) .yj. dm_eye(m, n-1, k, is_global)

    ! call dm_setvalues(L1, [(m-1)], [(n-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    ! res = L1 + L2
    ! call dm_destroy(L1, ierr)
    ! call dm_destroy(L2, ierr)

    call dm_set_implicit(res, ierr)
  end function OP_DZB


  ! % OP_DXF1_XY=[ -1  1  0  0  0  0  0]    
  ! %            [  0 -1  1  0  0  0  0]    
  ! %            [  0  0 -1  1  0  0  0]    
  ! %            [  0  0  0 -1  1  0  0]    
  ! %            [  0  0  0  0 -1  1  0]    
  ! %            [  0  0  0  0  0 -1  1]    
  ! %            [  0  0  0  0  0  0  0]
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

  
  ! % OP_DYF1_XY=[ -1  0  0  0  0 ]        
  ! %            [  1 -1  0  0  0 ]        
  ! %            [  0  1 -1  0  0 ]        ]
  ! %            [  0  0  1 -1  0 ]        
  ! %            [  0  0  0  1  0 ]        
  ! %                                      
  ! %                                      
  
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


  function OP_DZF(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    ! L1 = dm_eye(m, n, k, is_global)
    ! L2 = dm_zeros(m, 1, k, is_global) .yj. dm_eye(m, n-1, k, is_global)

    ! call dm_setvalues(L1, [(m-1)], [(n-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    ! res = L1 + L2
    ! call dm_destroy(L1, ierr)
    ! call dm_destroy(L2, ierr)

    call dm_set_implicit(res, ierr)
  end function OP_DZF


  function AXF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = MAT_AXF * A
    call dm_set_implicit(res, ierr)
  end function 

  function AXB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = MAT_AXB * A
    call dm_set_implicit(res, ierr)
  end function 

  
  function AYF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = A * MAT_AYF
    call dm_set_implicit(res, ierr)
  end function 

  function AYB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = A * MAT_AYB
    call dm_set_implicit(res, ierr)
  end function 

  function AZF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    ! res = A * MAT_AYB
    ! call dm_set_implicit(res, ierr)
  end function 

  function AZB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    ! res = A * MAT_AYB
    ! call dm_set_implicit(res, ierr)
  end function 
  
  function DXF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = MAT_DXF * A
    call dm_set_implicit(res, ierr)
  end function 

  function DXB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = MAT_DXB * A
    call dm_set_implicit(res, ierr)
  end function 

  function DXC(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = MAT_DXB * MAT_AXF * A
    call dm_set_implicit(res, ierr)
  end function 
  
  function DYF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = A * MAT_DYF
    call dm_set_implicit(res, ierr)
  end function 

  function DYB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = A * MAT_DYB
    call dm_set_implicit(res, ierr)
  end function 

  function DYC(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    res = A * MAT_AYF * MAT_DYB
    call dm_set_implicit(res, ierr)
  end function 
  
  function DZF(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    ! res = A * MAT_AYB
    ! call dm_set_implicit(res, ierr)
  end function 

  function DZB(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    ! res = A * MAT_AYB
    ! call dm_set_implicit(res, ierr)
  end function 

  function DZC(A) result(res)
    type(Matrix), intent(in) :: A
    type(Matrix) :: res
    integer :: ierr
    
    ! res = A * MAT_AYB
    ! call dm_set_implicit(res, ierr)
  end function 
  
  end module dm_op
