#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscdmdef.h> 

module ot_type

  type range
     integer :: lower = 0
     integer :: upper = 0
  end type range
  
  type box_info
     ! integer :: starts(3) = 0
     ! integer :: ends(3)   = 0
     type(range) :: rx
     type(range) :: ry
     type(range) :: rz     
  end type box_info
  
  type node_ptr
     type(node), pointer :: ptr => null()
  end type node_ptr

  type array_ptr
     type(array), pointer :: ptr => null()
  end type array_ptr

  !how the array distributed along x,y,z-directions
  type dist_info
     integer, allocatable, dimension(:) :: lx
     integer, allocatable, dimension(:) :: ly
     integer, allocatable, dimension(:) :: lz
     
     integer, allocatable, dimension(:) :: clx
     integer, allocatable, dimension(:) :: cly
     integer, allocatable, dimension(:) :: clz     
  end type dist_info
  
  type ref_info
     !0 for range, 1 for indices
     integer :: ref_index_type_x
     integer :: ref_index_type_y
     integer :: ref_index_type_z     

     !used only when ref_index_type = 0
     ! type(box_info) :: ref_box
     ! type(box_info) :: ref_lbox !local box

     type(range) :: range_x
     type(range) :: range_y
     type(range) :: range_z
     
     !used only when ref_index_type = 1
     integer, allocatable, dimension(:) :: iarr_x
     integer, allocatable, dimension(:) :: iarr_y
     integer, allocatable, dimension(:) :: iarr_z     

     !0 for array, 1 for node
     integer :: ref_target_type
     ! type(node),  pointer :: ref_node   => null()
     ! type(array),pointer :: ref_array => null()
  end type ref_info
  
  type node
     type(array), pointer :: data => null()

     !the operands for a operation, e.g. A + B,, A and B are operands
     !for a binary operator, the number operands is always 2
     type(node_ptr), allocatable, dimension(:) :: operands

     !node type, 0 for data, 1 for '+', 2 for '-',
     !3 for '*', 4 for '/', etc.
     integer :: node_type

     !only used when node_type is 'type_ref'
     type(ref_info) :: ref

     ! the node shape inherit from its left/right
     ! expression or the array within it
     integer :: m_shape(3) = 0

     real(8) :: alpha = 0
     real(8) :: beta  = 1
     real(8) :: scalar = 0

     !number of arguments
     integer :: num_args = 0
     
     !arguments for the node
     real(8) :: args(8) = 0

     !reference counter
     integer :: ref_cnt = 0

     !default is lvalue
     character(len=1) :: var_type = 'l'

     !store its hash
     integer(kind=8) :: hash = 0

     !the unique id for each node
     integer :: id
  end type node

  type array
     !dim and shape should not be specified by user!
     !integer :: m_dim = 0
     integer :: m_shape(3) = (/0,0,0/)

     Vec :: data = 0

     logical :: is_field = .false.
     
     integer :: grid_pos = -1

     !local data block
     type(box_info) :: local_block

     !reference counter
     integer :: ref_cnt = 0

     !default is lvalue, lvalue, rvalue
     character(len=1) :: var_type = 'l'
  end type array

  type grid
     character(1) :: grid_type !'B' or 'C'
     type(array), pointer :: dx_3d => null()
     type(array), pointer :: dy_3d => null()
     type(array), pointer :: dz_3d => null()     
  end type grid

  
end module
