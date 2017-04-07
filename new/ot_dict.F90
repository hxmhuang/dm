#include "type.h"

module ot_dict
  use ot_common
  
  type dict_item
     integer(kind=8) :: key = 0
     C_POINTER :: func = 0
     type(dict_item), pointer :: left => null()
     type(dict_item), pointer :: right => null()
  end type dict_item

  interface disp
     module procedure dict_disp
  end interface disp
contains

  recursive subroutine dict_disp(dict)
    implicit none
    type(dict_item), intent(in) :: dict

    write(*,'(I0.1, A, Z16.16)'), dict%key, " : 0x", dict%func
    if(associated(dict%left))  call dict_disp(dict%left)
    if(associated(dict%right)) call dict_disp(dict%right)    
  end subroutine
  
  recursive subroutine dict_add(dict, key, func)
    implicit none
    integer(kind=8) :: key
    C_POINTER :: func
    type(dict_item), pointer, intent(inout) :: dict
    type(dict_item), pointer :: left, right
    
    if(.not. associated(dict)) allocate(dict)

    call assert(key /= dict%key, __FILE__, __LINE__,&
         'Error: key already exists.')
    
    if(key < dict%key) then
       if(associated(dict%left)) then
          call dict_add(dict%left, key, func)
       else
          allocate(dict%left)
          dict%left%key = key
          dict%left%func = func
       endif
    else
       if(associated(dict%right)) then
          call dict_add(dict%right, key, func)
       else
          allocate(dict%right)
          dict%right%key = key
          dict%right%func = func
       end if
    end if
  end subroutine

  recursive function dict_get(dict, key) result(func)
    implicit none
    integer(kind=8) :: key
    C_POINTER :: func
    type(dict_item), pointer, intent(in) :: dict

    if(.not. associated(dict)) then
       func = 0
       return
    end if
    
    if(dict%key == key) then
       func = dict%func
       return
    endif

    if(key <= dict%key) then
       func = dict_get(dict%left, key)
    else
       func = dict_get(dict%right, key)
    end if
  end function

end module ot_dict
