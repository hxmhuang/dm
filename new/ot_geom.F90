
module ot_geom
  use ot_type
  
  interface operator (.eq.)
     module procedure box_eq
     module procedure range_eq
  end interface operator (.eq.)
  
  interface operator (.and.)
     module procedure box_and
  end interface operator (.and.)

  interface operator (.in.)
     module procedure box_in
  end interface operator (.in.)
  
  interface shape
     module procedure box_shape
  end interface shape

  interface disp
     module procedure disp_box
  end interface disp

  type(range) :: range_all
contains
  subroutine init_geom()
    implicit none
    range_all%lower = -huge(int(0))
    range_all%lower = huge(int(0))
  end subroutine

  function range_eq(a, b) result(res)
    implicit none
    type(range), intent(in) :: a
    type(range), intent(in) :: b
    logical :: res
        
    res = (a%lower == b%lower) .and. &
         (a%upper == b%upper)
  end function

  subroutine range_to_box(b, rx, ry, rz)
    implicit none
    type(box_info), intent(out) :: b
    type(range), intent(in) :: rx, ry, rz
    
    b%starts(1) = rx%lower
    b%starts(2) = ry%lower
    b%starts(3) = rz%lower
    b%ends(1) = rx%upper
    b%ends(2) = ry%upper
    b%ends(3) = rz%upper
  end subroutine

  subroutine range_to_ref(b, rx, ry, rz)
    implicit none
    type(ref_info), intent(out) :: b
    type(range), intent(in) :: rx, ry, rz
    
    b%range_x = rx
    b%range_y = ry
    b%range_z = rz

    b%ref_index_type_x = 0
    b%ref_index_type_y = 0
    b%ref_index_type_z = 0
  end subroutine
  
  subroutine disp_box(o, prefix)
    implicit none
    type(box_info) :: o
    character(len=*),optional :: prefix
    
    if(present(prefix)) then
       write(*, "(2X, A)") prefix
    else
       write(*, "(2X, A)") "ans =  "
    endif
    
    write(*, "(6X, 3I4.1)") o%starts
    write(*, "(6X, 3I4.1)") o%ends
  end subroutine
  
  function r(a, b) result(res)
    implicit none
    integer :: a
    integer :: b
    type(range) :: res
    res%lower = a
    res%upper = b
  end function

  !> if box a is inside box b
  function box_in(a, b) result(res)
    implicit none
    type(box_info), intent(in) :: a, b
    logical :: res

    res = (.not. any(a%starts < b%starts)) &
         .and. (.not. any(a%ends > b%ends))
  end function
  
  function has_intersection(a, b) result(res)
    implicit none
    type(box_info), intent(in) :: a, b
    logical :: res
    
    res = .not. any(shape(a .and. b) < 1)
    
  end function

  function is_valid(a) result(res)
    implicit none
    type(box_info), intent(in) :: a
    logical :: res

    res = any(shape(a) <= 0)
    
  end function
  
  function box_shape(a) result(res)
    implicit none
    type(box_info) :: a
    integer, dimension(3) :: res
    res = a%ends - a%starts + 1
  end function
  
  function box_and(a, b) result(res)
    implicit none
    type(box_info), intent(in) :: a, b
    type(box_info) :: res

    res%starts(1) = max(a%starts(1), b%starts(1))
    res%ends(1) = min(a%ends(1), b%ends(1))

    res%starts(2) = max(a%starts(2), b%starts(2))
    res%ends(2) = min(a%ends(2), b%ends(2))

    res%starts(3) = max(a%starts(3), b%starts(3))
    res%ends(3) = min(a%ends(3), b%ends(3))
    
  end function
  
  function box_eq(a, b) result(res)
    implicit none
    type(box_info), intent(in) :: a
    type(box_info), intent(in) :: b
    logical :: res

    res = (a%starts(1) == b%starts(1)) .and. &
         (a%starts(2) == b%starts(2)) .and. &
         (a%starts(3) == b%starts(3)) .and. &
         (a%ends(1) == b%ends(1)) .and. &
         (a%ends(2) == b%ends(2)) .and. &
         (a%ends(3) == b%ends(3))
  end function
  
end module
