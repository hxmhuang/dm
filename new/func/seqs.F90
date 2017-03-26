module mod_seqs
  use ot_tensor
  use ot_data
  
  interface seqs
     module procedure seqs_1d
     module procedure seqs_2d
     module procedure seqs_3d
  end interface seqs
  
contains
  function seqs_1d(m) result(res)
    integer :: m
    type(tensor) :: res
    integer :: ierr
    
    call tensor_set_shape(res, (/m/))
    call data_seqs(res%data, res%m_shape, ierr)
  end function

  function seqs_2d(m, n) result(res)
    integer :: m, n
    type(tensor) :: res

    call tensor_set_shape(res, (/m,n/))
    call data_seqs(res%data, res%m_shape, ierr)
  end function

  function seqs_3d(m, n, k) result(res)
    integer :: m, n, k
    type(tensor) :: res

    call tensor_set_shape(res, (/m,n,k/))
    call data_seqs(res%data, res%m_shape, ierr)
  end function
end module

