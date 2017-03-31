
module ot_geom
  use ot_type
  
  interface operator (.eq.)
     module procedure box_eq
     module procedure range_eq
     module procedure dist_eq
  end interface operator (.eq.)

  interface operator (.ne.)
     module procedure dist_ne
  end interface operator (.ne.)
  
  interface operator (.and.)
     module procedure box_and
  end interface operator (.and.)

  interface operator (.in.)
     module procedure box_in
     module procedure range_in
  end interface operator (.in.)
  
  interface shape
     module procedure box_shape
  end interface shape

  interface disp
     module procedure disp_box
  end interface disp

  interface operator(-)
     module procedure range_minus_int
  end interface operator(-)

  interface operator(+)
     module procedure range_plus_int
  end interface operator(+)

  type(range) :: range_all
contains
  subroutine init_geom(ierr)
    implicit none
    integer, intent(out) :: ierr
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

  subroutine box_to_indices(idx, box, gshape)
    implicit none
    integer, allocatable, intent(out) :: idx(:)
    type(box_info) :: box, box1
    integer :: gshape(3) !global shape
    integer :: gm, gn, bs(3), i, k, j, cnt

    if(.not. is_valid(box)) then
       allocate(idx(0))
       return
    end if
    
    bs = shape(box)
    allocate(idx(bs(1) * bs(2) * bs(3)))
    
    gm = gshape(1)
    gn = gshape(2)
    
    cnt = 1

    do k = box%rz%lower, box%rz%upper
       do j = box%ry%lower, box%ry%upper
          do i = box%rx%lower, box%rx%upper
             idx(cnt) = i + j * gm + k * gm * gn
             cnt = cnt + 1
          end do
       end do
    end do
  end subroutine
  
  subroutine box_display(box, prefix)
    implicit none
    type(box_info) :: box
    character(len=*), optional :: prefix

    if(present(prefix)) then
       write(*, "(2X,A)") prefix
    else
       write(*, "(2X,A)") "ans = "
    end if

    write(*, "(4X, 3g10.5)") box%rx%lower, &
         box%ry%lower,box%rz%lower
    write(*, "(4X, 3g10.5)") box%rx%upper, &
         box%ry%upper,box%rz%upper
  end subroutine
  
  function dist_ne(a, b) result(res)
    implicit none
    type(dist_info), intent(in) :: a
    type(dist_info), intent(in) :: b
    logical :: res
    
    res = .not. (a .eq. b)
  end function
  
  function dist_eq(a, b) result(res)
    implicit none
    type(dist_info), intent(in) :: a
    type(dist_info), intent(in) :: b
    logical :: res
    integer :: nxa, nya, nza, nxb, nyb, nzb, i

    nxa = size(a%lx)
    nya = size(a%ly)
    nza = size(a%lz)
    nxb = size(b%lx)
    nyb = size(b%ly)
    nzb = size(b%lz)

    !check the dimension first
    if(nxa /= nxb .or. nya /= nyb .or. nza /= nzb) then
       res = .false.
       return
    end if

    do i = 1, nxa
       if(a%lx(i) /= b%lx(i)) then
          res = .false.
          return
       end if
    enddo
    
    do i = 1, nya
       if(a%ly(i) /= b%ly(i)) then
          res = .false.
          return
       end if
    enddo
    
    do i = 1, nza
       if(a%lz(i) /= b%lz(i)) then
          res = .false.
          return
       end if
    enddo
    
    res = .true.
  end function
  
  subroutine range_to_box(b, rx, ry, rz)
    implicit none
    type(box_info), intent(out) :: b
    type(range), intent(in) :: rx, ry, rz
    
    b%rx%lower = rx%lower
    b%ry%lower = ry%lower
    b%rz%lower = rz%lower
    b%rx%upper = rx%upper
    b%ry%upper = ry%upper
    b%rz%upper = rz%upper
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

  function shape_to_ref(s) result (res)
    implicit none
    integer, intent(in) :: s(3)
    type(ref_info)  :: res
    
    res%range_x = r(0, s(1)-1)
    res%range_y = r(0, s(2)-1)
    res%range_z = r(0, s(3)-1)
    
    res%ref_index_type_x = 0
    res%ref_index_type_y = 0
    res%ref_index_type_z = 0    
  end function
  
  function range_minus_int(rgn, num) result(res)
    implicit none
    type(range), intent(in) :: rgn
    integer, intent(in) :: num
    type(range) :: res

    res%upper = rgn%upper - num
    res%lower = rgn%lower - num    
  end function

  function range_plus_int(rgn, num) result(res)
    implicit none
    type(range), intent(in) :: rgn
    integer, intent(in) :: num
    type(range) :: res

    res%upper = rgn%upper + num
    res%lower = rgn%lower + num    
  end function
  
  subroutine disp_box(o, prefix)
    implicit none
    type(box_info) :: o
    character(len=*),optional :: prefix
    
    if(present(prefix)) then
       write(*, "(2X, A)") prefix
    else
       write(*, "(2X, A)") "ans =  "
    endif
    
    write(*, "(6X, 3I4.1)") o%rx%lower, &
         o%ry%lower,o%rz%lower
    write(*, "(6X, 3I4.1)") o%rx%upper, &
         o%ry%upper,o%rz%upper        
  end subroutine
  
  function r(a, b) result(res)
    implicit none
    integer :: a
    integer :: b
    type(range) :: res
    res%lower = a
    res%upper = b
  end function

  !> check if range a is inside range b
  function range_in(a, b) result(res)
    implicit none
    type(range), intent(in) :: a, b
    logical res
    res = (a%upper <= b%upper) .and. (a%lower >= b%lower)
  end function
  
  !> if box a is inside box b
  function box_in(a, b) result(res)
    implicit none
    type(box_info), intent(in) :: a, b
    logical :: res

    res = (a%rx .in. b%rx) .and. &
         (a%ry .in. b%ry) .and. (a%rz .in. b%rz)
    
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

    res = .not. any(shape(a) <= 0)
    
  end function

  !> the shape of box
  function box_shape(a) result(res)
    implicit none
    type(box_info) :: a
    integer, dimension(3) :: res
    !res = a%ends - a%starts + 1
    res(1) = a%rx%upper - a%rx%lower + 1
    res(2) = a%ry%upper - a%ry%lower + 1
    res(3) = a%rz%upper - a%rz%lower + 1    
  end function

  !> calculate the cross section of two boxes
  function box_and(a, b) result(res)
    implicit none
    type(box_info), intent(in) :: a, b
    type(box_info) :: res

    res%rx%lower = max(a%rx%lower, b%rx%lower)
    res%rx%upper = min(a%rx%upper, b%rx%upper)

    res%ry%lower = max(a%ry%lower, b%ry%lower)
    res%ry%upper = min(a%ry%upper, b%ry%upper)

    res%rz%lower = max(a%rz%lower, b%rz%lower)
    res%rz%upper = min(a%rz%upper, b%rz%upper)
  end function

  subroutine get_ref_box(sub_box2, box2, sub_box1, box1)
    implicit none
    type(box_info), intent(in) :: box1, box2, sub_box1
    type(box_info), intent(out) :: sub_box2

    sub_box2%rx = sub_box1%rx - box1%rx%lower + box2%rx%lower
    sub_box2%ry = sub_box1%ry - box1%ry%lower + box2%ry%lower
    sub_box2%rz = sub_box1%rz - box1%rz%lower + box2%rz%lower    

  end subroutine
  
  function box_eq(a, b) result(res)
    implicit none
    type(box_info), intent(in) :: a
    type(box_info), intent(in) :: b
    logical :: res

    res = (a%rx%lower == b%rx%lower) .and. &
         (a%ry%lower == b%ry%lower) .and. &
         (a%rz%lower == b%rz%lower) .and. &
         (a%rx%upper == b%rx%upper) .and. &
         (a%ry%upper == b%ry%upper) .and. &
         (a%rz%upper == b%rz%upper)
  end function
  
end module
