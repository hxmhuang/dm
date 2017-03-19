
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 
#:include "type_def.fypp"
#include "node_type.h"

module ot_ref
  use ot_geom
  use petsc_helper
  ! use ot_tensor
  ! use ot_node
  ! type tensor_ptr_array
  !    integer :: data_size = 0
  !    type(tensor_ptr), allocatable, dimension(:) :: v
  ! end type tensor_ptr_array

  

  ! interface size
  !    module procedure size_arr
  ! end interface size

  !private :: shape 

contains
  ! function size_arr(arr) result(res)
  !   implicit none
  !   type(tensor_ptr_array), intent(in) :: arr
  !   integer :: res

  !   res = arr%data_size
  ! end function

  ! subroutine test_size()
  !   implicit none
  !   type(tensor_ptr_array) :: p
  !   integer :: a(102)
  !   allocate(p%v(100))
  !   p%data_size = 100    
  !   print*, "size=", size(p)
  !   print*, "size_a = ", size(a)
  ! end subroutine



end module
