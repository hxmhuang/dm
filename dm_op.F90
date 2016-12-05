module dm_op
  use dm_type
  use dm
  implicit none
#include "mat_type.h"
contains

  function OP_AXF(m, n, k, isGlobal) result(res)
    integer, intent(in) :: m, n, k
    logical, intent(in), optional :: isGlobal
    type(Matrix) :: res
    type(Matrix) :: L1, L2
    logical :: is_global = .true.
    integer :: ierr
    integer :: i

    if(present(isGlobal)) is_global = isGlobal
    L1 = dm_eye(m, n, k, is_global)
    L2 = dm_zeros(m, 1, k, is_global) .yj. dm_eye(m, n-1, k, is_global)

    call dm_setvalues(L1, [(m-1)], [(n-1)], [(i, i=0,k-1)], [(2,i=0,k-1)], ierr)

    res = L1 + L2
    call dm_destroy(L1, ierr)
    call dm_destroy(L2, ierr)

    call dm_set_implicit(res, ierr)
  end function OP_AXF

  function OP_AYF(m, n, k, isGlobal) result(res)
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
  end function OP_AYF

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


  function OP_AXB(m, n, k, isGlobal) result(res)
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
  end function OP_AXB


  function OP_AYB(m, n, k, isGlobal) result(res)
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


  function OP_DXB(m, n, k, isGlobal) result(res)
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
  end function OP_DXB


  function OP_DYB(m, n, k, isGlobal) result(res)
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


  function OP_DXF(m, n, k, isGlobal) result(res)
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
  end function OP_DXF


  function OP_DYF(m, n, k, isGlobal) result(res)
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


  function OP_DXC(m, n, k, isGlobal) result(res)
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
  end function OP_DXC


  function OP_DYC(m, n, k, isGlobal) result(res)
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
  end function OP_DYC


  function OP_DZC(m, n, k, isGlobal) result(res)
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
  end function OP_DZC

end module dm_op
