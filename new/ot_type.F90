
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 
#:include "type_def.fypp"
#include "node_type.h"

module ot_type
  use ot_geom
  
  type grid
     character(1) :: grid_type !'B' or 'C'
     type(tensor), pointer :: dx_3d => null()
     type(tensor), pointer :: dy_3d => null()
     type(tensor), pointer :: dz_3d => null()     
  end type grid

  type tensor
     !dim and shape should not be specified by user!
     integer :: m_dim = 0, m_shape(3) = (/0,0,0/)

     Vec :: data = 0

     !check if it is an implicit variable
     logical :: is_implicit = .false.

     logical :: is_field = .false.
     
     integer :: grid_pos

     type(block) :: local_block
     
  end type tensor

  type tensor_ptr
     type(tensor), pointer :: ptr => null()
  end type tensor_ptr

  type tensor_ptr_array
     integer :: data_size = 0
     type(tensor_ptr), allocatable, dimension(:) :: v
  end type tensor_ptr_array
  
  type node_ptr
     type(node), pointer :: ptr => null()
  end type node_ptr
  
  type node
     type(tensor), pointer :: data => null()

     type(block) :: local_block
     
     !the operands for a operation, e.g. A + B,, A and B are operands
     !for a binary operator, the number operands is always 2
     type(node_ptr), allocatable, dimension(:) :: operands

     !optimized operands, e.g. A + B + C, opt_operands are (A,B,C)
     type(node_ptr), allocatable, dimension(:) :: opt_operands

     !node type, 0 for data, 1 for '+', 2 for '-',
     !3 for '*', 4 for '/', etc.
     integer :: node_type

     ! the node shape inherit from its left/right
     ! expression or the tensor within it
     integer :: m_dim = 0
     integer :: m_shape(3) = (/0,0,0/)

     real(8) :: alpha = 0
     real(8) :: beta  = 1
     real(8) :: scalar = 0

     logical :: is_implicit

     !the unique id for each node
     integer :: id

     !arguments for the node
     real(8) :: args(10)
  end type node

  interface size
     module procedure size_arr
  end interface size
  
contains
  function size_arr(arr) result(res)
    implicit none
    type(tensor_ptr_array), intent(in) :: arr
    integer :: res

    res = arr%data_size
  end function

  subroutine test_size()
    implicit none
    type(tensor_ptr_array) :: p
    allocate(p%v(100))
    print*, "size=", size(p)
  end subroutine
    
end module
