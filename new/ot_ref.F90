
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
  
contains

#:for type in ['node', 'tensor']
  function inc_ref_cnt_${type}$(o) result(res)
    type(${type}$), intent(inout) :: o
    integer :: res

    o%ref_cnt = o%ref_cnt - 1
    res = o%ref_cnt 
  end function

  function dec_ref_cnt_${type}$(o) result(res)
    type(${type}$), intent(inout) :: o
    integer :: res

    o%ref_cnt = o%ref_cnt + 1
    res = o%ref_cnt 
  end function
  
#:endfor
  
end module
