module mod_ones
  use ot_tensor
  use ot_data
  
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
    
    call tensor_set_shape(res, (/m/))
    call data_consts(res%data, real(1, 8), res%m_shape, ierr)
  end function

  function ones_2d(m, n) result(res)
    integer :: m, n
    type(tensor) :: res
    real(kind=8) :: one = 1

    call tensor_set_shape(res, (/m,n/))        
    call data_consts(res%data, real(1, 8), res%m_shape, ierr)
  end function

  function ones_3d(m, n, k) result(res)
    integer :: m, n, k
    type(tensor) :: res
    real(kind=8) :: one = 1

    call tensor_set_shape(res, (/m,n,k/))        
    call data_consts(res%data, real(1, 8), res%m_shape, ierr)
  end function
end module

