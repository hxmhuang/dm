
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 
#:include "type_def.fypp"
#include "type.h"

module ot_ref
  use ot_type

  interface inc_ref_cnt
#:for type in ['node', 'tensor']
     module procedure inc_ref_cnt_${type}$
#:endfor     
  end interface inc_ref_cnt
  
  interface dec_ref_cnt
#:for type in ['node', 'tensor']
     module procedure dec_ref_cnt_${type}$
#:endfor     
  end interface dec_ref_cnt

  interface is_rvalue
     module procedure is_rvalue_node
     module procedure is_rvalue_tensor     
  end interface is_rvalue
  
  interface is_lvalue
     module procedure is_lvalue_node
     module procedure is_lvalue_tensor     
  end interface is_lvalue

  interface set_rvalue
     module procedure set_rvalue_node
     module procedure set_rvalue_tensor
  end interface set_rvalue

  interface set_lvalue
     module procedure set_lvalue_node
     module procedure set_lvalue_tensor
  end interface set_lvalue
  
contains

#:for type in ['node', 'tensor']
  function inc_ref_cnt_${type}$(o) result(res)
    type(${type}$), intent(inout) :: o
    integer :: res

    o%ref_cnt = o%ref_cnt + 1
    res = o%ref_cnt 
  end function

  function dec_ref_cnt_${type}$(o) result(res)
    type(${type}$), intent(inout) :: o
    integer :: res

    o%ref_cnt = o%ref_cnt - 1
    res = o%ref_cnt 
  end function
  
#:endfor

#:for type in ['node', 'tensor']
  function is_rvalue_${type}$(o) result(res)
    implicit none
    type(${type}$), intent(in) :: o
    logical :: res
    
    res = (o%var_type == 'r')
  end function

  function is_lvalue_${type}$(o) result(res)
    implicit none
    type(${type}$), intent(in) :: o
    logical :: res
    
    res = (o%var_type == 'l')
  end function

  subroutine set_rvalue_${type}$(o)
    implicit none
    type(${type}$), intent(inout) :: o
    
    o%var_type = 'r'
  end subroutine

  subroutine set_lvalue_${type}$(o)
    implicit none
    type(${type}$), intent(inout) :: o
    
    o%var_type = 'l'
  end subroutine
#:endfor
  
end module
