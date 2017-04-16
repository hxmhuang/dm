module ot_func
  use ot_tensor
  use ot_data
  
  interface consts
     module procedure consts_int
     module procedure consts_real
     module procedure consts_real8
  end interface consts
  
contains

#:for type in [['int', 'integer'], ['real', 'real'],['real8','real(8)']]
  function consts_${type[0]}$(val, m, opt_n, opt_k) result(res)
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
    call set_rvalue(res)
  end function
  
#:endfor

  function ones(m, opt_n, opt_k) result(res)
    implicit none
#include "petsc.h"
    integer, intent(in) :: m
    integer, intent(in), optional :: opt_n, opt_k
    type(tensor), pointer :: res
    integer :: n, k, ierr
    Vec :: data
    
    allocate(res)
    
    n = 1; k = 1;
    
    if(present(opt_n)) n = opt_n
    if(present(opt_k)) k = opt_k

    !call tensor_set_shape(res, (/m,n,k/))
    call data_consts(data, real(1,8), (/m, n, k/), ierr)
    call tensor_bind_data(res, data)
    call set_rvalue(res)
  end function

  function seqs(m, opt_n, opt_k) result(res)
    implicit none
#include "petsc.h"
    integer, intent(in) :: m
    integer, intent(in), optional :: opt_n, opt_k
    type(tensor), pointer :: res
    integer :: n, k, ierr
    Vec :: data
    allocate(res)
    
    n = 1; k = 1;
    
    if(present(opt_n)) n = opt_n
    if(present(opt_k)) k = opt_k

    !call tensor_set_shape(res, (/m,n,k/))
    !call data_seqs(res%data, res%m_shape, ierr)
    call data_seqs(data, (/m,n,k/), ierr)
    call tensor_bind_data(res, data)
    call set_rvalue(res)
  end function

  function rand(m, opt_n, opt_k) result(res)
    implicit none
#include "petsc.h"
    integer, intent(in) :: m
    integer, intent(in), optional :: opt_n, opt_k
    type(tensor), pointer :: res
    integer :: n, k, ierr
    Vec :: data
    allocate(res)
    
    n = 1; k = 1;
    
    if(present(opt_n)) n = opt_n
    if(present(opt_k)) k = opt_k

    !call tensor_set_shape(res, (/m,n,k/))
    ! call data_rand(res%data, res%m_shape, ierr)
    call data_rand(data, (/m,n,k/), ierr)
    call tensor_bind_data(res, data)
    call set_rvalue(res)
  end function
  
end module

