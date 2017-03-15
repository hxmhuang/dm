
module ot_geom
  type block
     integer :: starts(3) = (/0,0,0/)
     integer :: counts(3) = (/0,0,0/)
  end type block

  interface operator(==)
     module procedure block_eq
  end interface operator(==)
  
contains
  function block_eq(a, b) result(res)
    implicit none
    type(block), intent(in) :: a
    type(block), intent(in) :: b
    logical :: res

    res = (a%starts(1) == b%starts(1)) .and. &
         (a%starts(2) == b%starts(2)) .and. &
         (a%starts(3) == b%starts(3)) .and. &
         (a%counts(1) == b%counts(1)) .and. &
         (a%counts(2) == b%counts(2)) .and. &
         (a%counts(3) == b%counts(3))
  end function
  
end module
