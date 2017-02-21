module dm_func_consts
  use dm_tensor
  use dm_data
  
  interface consts
     module procedure consts_1d
     module procedure consts_2d
     module procedure consts_3d
  end interface consts
  
contains
  function consts_1d(val,m) result(res)
    integer :: m
    real(kind=8) :: val
    type(tensor) :: res
    integer :: ierr
    
    res = tensor_new((/m/))
    call data_constants(res%data, val, (/m/), ierr)
  end function

  function consts_2d(val,m, n) result(res)
    integer :: m, n
    real(kind=8) :: val    
    type(tensor) :: res
    print*, "val=", val
    res = tensor_new((/m, n/))
    call data_constants(res%data,  val, (/m,n/), ierr)
  end function

  function consts_3d(val,m, n, k) result(res)
    integer :: m, n, k
    real(kind=8) :: val    
    type(tensor) :: res
    
    res = tensor_new((/m,n,k/))
    call data_constants(res%data,  val, (/m,n,k/), ierr)
  end function
end module

