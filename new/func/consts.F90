module mod_consts

  use ot_tensor
  use ot_data
  
  interface consts
     module procedure consts_3d_int
     module procedure consts_3d_real
     module procedure consts_3d_real8
  end interface consts
  
contains

#:for type in [['int', 'integer'], ['real', 'real'],['real8','real(8)']]
  function consts_3d_${type[0]}$(val, m, opt_n, opt_k) result(res)
    implicit none
    ${type[1]}$ :: val
    integer, intent(in) :: m
    integer, intent(in), optional :: opt_n, opt_k
    type(tensor) :: res
    integer :: n, k, ierr

    n = 1; k = 1;
    
    if(present(opt_n)) n = opt_n
    if(present(opt_k)) k = opt_k

    call tensor_set_shape(res, (/m,n,k/))
    call data_consts(res%data, real(val,8), res%m_shape, ierr)    
    res%var_type = 'r'    
  end function
  
#:endfor
end module

