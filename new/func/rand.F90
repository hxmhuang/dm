module dm_func_rand
  use ot_tensor
  use ot_data
  
  interface rand
     module procedure rand_3d
  end interface rand
  
contains

  function rand_3d(m, opt_n, opt_k) result(res)
    implicit none
    integer, intent(in) :: m
    integer, intent(in), optional :: opt_n, opt_k
    type(tensor) :: res
    integer :: n, k, ierr

    n = 1; k = 1;
    
    if(present(opt_n)) n = opt_n
    if(present(opt_k)) k = opt_k

    call tensor_set_shape(res, (/m,n,k/))
    call data_rand(res%data, res%m_shape, ierr)
    res%var_type = 'r'    
  end function

end module

