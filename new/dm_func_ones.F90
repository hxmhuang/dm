module dm_func_ones
  use dm_tensor
  use dm_data
  
  interface ones
     module procedure ones_1d
     module procedure ones_2d
     module procedure ones_3d
  end interface ones
  
contains
  function ones_1d(m) result(res)
    integer :: m
    type(tensor) :: res
    real(kind=8) :: one = 1
    integer :: ierr
    
    res = tensor_new((/m/))
    call data_constants(res%data, one, (/m/), ierr)
  end function

  function ones_2d(m, n) result(res)
    integer :: m, n
    type(tensor) :: res
    real(kind=8) :: one = 1
    
    res = tensor_new((/m, n/))
    call data_constants(res%data,  one, (/m,n/), ierr)
  end function

  function ones_3d(m, n, k) result(res)
    integer :: m, n, k
    type(tensor) :: res
    real(kind=8) :: one = 1
    
    res = tensor_new((/m,n,k/))
    call data_constants(res%data,  one, (/m,n,k/), ierr)
  end function
end module

