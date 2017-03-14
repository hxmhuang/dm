module dm_func_set
  use ot_tensor
  use ot_data
  
  interface set
     module procedure set1
     module procedure set2     
  end interface set
  
contains
  function set(dst, val) result(res)
    type(tensor), intent(inout) :: dst
    real*8 :: val
    real(kind=8) :: one = 1
    integer :: ierr

    call data_set(dst, val)
    call data_constants(res%data, one, (/m/), ierr)
  end function

  function set(dst, src) result(res)
    type(tensor), intent(inout) :: dst
    real*8 :: val
    real(kind=8) :: one = 1
    integer :: ierr

    call data_set(dst, val)
    call data_constants(res%data, one, (/m/), ierr)
  end function

end module

