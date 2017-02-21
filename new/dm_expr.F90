module dm_expr
  use dm_common
  use dm_tensor
#include "dm_type.h"
  
  public :: operator(+),operator(-),operator(*),operator(/)
  public :: assignment(=)

  interface operator (+)
     module procedure node_plus
  end interface 

  interface operator (-)
     module procedure node_minus
  end interface 

  interface operator (*)
     module procedure node_mult
  end interface 

  interface operator (/)
     module procedure node_divd
  end interface 
  
  interface assignment (=)
     module procedure tensor_assign
  end interface assignment (=)

contains
  
  function node_plus(A, B) result(C)
    implicit none
    type(tensor), intent(in), target :: A
    type(tensor), intent(in), target :: B
    type(tensor) :: C

    C = create_tensor_node(type_plus, A, B);        
  end function node_plus

  function node_minus(A, B) result(C)
    implicit none
    type(tensor), intent(in), target :: A
    type(tensor), intent(in), target :: B
    type(tensor) :: C

    C = create_tensor_node(type_minus, A, B);    
  end function 
  
  function node_mult(A, B) result(C)
    implicit none
    type(tensor), intent(in), target :: A
    type(tensor), intent(in), target :: B
    type(tensor) :: C

    C = create_tensor_node(type_mult, A, B);    
  end function 

  function node_divd(A, B) result(C)
    implicit none
    type(tensor), intent(in), target :: A
    type(tensor), intent(in), target :: B
    type(tensor) :: C

    C = create_tensor_node(type_divd, A, B);
  end function 

  
  function create_tensor_node(node_type, left, right) result(A)
    implicit none

    type(tensor), intent(in), target :: left, right
    integer :: node_type
    type(tensor) :: A

    call assert(left%m_dim == right%m_dim, __FILE__,  __LINE__)
    
    A%m_dim   = left%m_dim
    A%m_shape = left%m_shape

    A%left => left
    A%right => right
    A%node_type = node_type
    A%is_implicit = .true.
  end function
  
  recursive subroutine eval(B)
    implicit none
    type(tensor), intent(inout) :: B

    if(B%node_type == type_data) return 

    call assert(associated(B%left) .and. associated(B%right), &
         __FILE__, __LINE__)
    
    !evaluate left and right leaves recursively.
    if(B%left%node_type  /= type_data) call eval(B%left)
    if(B%right%node_type /= type_data) call eval(B%right)       

    select case(B%node_type)
    case (type_plus)
       call tensor_plus(B, B%left, B%right)  !C = A + B
    case (type_minus)
       call tensor_minus(B, B%left, B%right) !C = A - B
    case (type_mult)
       call tensor_mult(B, B%left, B%right)  !C = A .* B
    case (type_divd)
       call tensor_divd(B, B%left, B%right)  !C = A ./ B
    end select
    
    B%node_type = type_data
  end subroutine eval
  
  subroutine tensor_assign(A, B)
    implicit none
    type(tensor), intent(out) :: A
    type(tensor), intent(in) :: B
    type(tensor) :: C

    
    call tensor_copy(A, B)
    if(B%node_type /= type_data) call eval(A)
    
    ! if(B%node_type == type_data) then
    !    call tensor_copy(A, B)
    ! else
    !    call tensor_copy(C, B)
    !    call eval(C)
    !    call tensor_copy(A, C)
    ! endif
    
  end subroutine tensor_assign


end module

