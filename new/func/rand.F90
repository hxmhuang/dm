module dm_func_rand
  use ot_tensor
  use ot_data
  
  interface rand
     module procedure rand_1d
     module procedure rand_2d
     module procedure rand_3d
  end interface rand
  
contains
  function rand_1d(m) result(res)
    integer :: m
    type(tensor) :: res
    real(kind=8) :: one = 1
    integer :: ierr
    
    call tensor_set_shape(res, (/m/))
    call data_rand(res%data, res%m_shape, ierr)
  end function

  function rand_2d(m, n) result(res)
    integer :: m, n
    type(tensor) :: res
    real(kind=8) :: one = 1

    call tensor_set_shape(res, (/m,n/))
    call data_rand(res%data, res%m_shape, ierr)
  end function

  function rand_3d(m, n, k) result(res)
    integer :: m, n, k
    type(tensor) :: res
    real(kind=8) :: one = 1

    call tensor_set_shape(res, (/m,n,k/))
    call data_rand(res%data, res%m_shape, ierr)
  end function
end module

