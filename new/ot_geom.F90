
module ot_geom
  use ot_type
  use ot_common
  use ot_buffer
  
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

  interface num_data
     module procedure box_num_data
  end interface num_data
  
  interface disp
     module procedure disp_box
     module procedure disp_ref
     module procedure disp_dist
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

  subroutine box_to_indices1(idx, box, gshape)
    implicit none
    integer, allocatable, intent(out) :: idx(:)
    type(box_info) :: box, box1
    integer :: gshape(3) !global shape
    integer :: gm, gn, bs(3), i, k, j, cnt
    integer :: ierr

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
  
  subroutine disp_ref(ref, prefix, indent)
    implicit none
    type(ref_info), intent(in) :: ref
    character(len=*), optional :: prefix
    character(len=*), optional :: indent

    if(present(prefix)) then
       write(*, "(4X, A)") prefix
    else
       write(*, "(4X, A)") "ref box : "
    endif

    if(ref%ref_index_type_x == 0) then
       write(*, "(10X, A, I0.1, A, I0.1, A)") &
            'm : [', ref%range_x%lower+1, ' : ', &
            ref%range_x%upper+1, ']'
    end if
    
    if(ref%ref_index_type_y == 0) then
       write(*, "(10X, A, I0.1, A, I0.1, A)") &
            'n : [', ref%range_y%lower+1, ' : ', &
            ref%range_y%upper+1, ']'
    end if
    
    if(ref%ref_index_type_y == 0) then
       write(*, "(10X, A, I0.1, A, I0.1, A)") &
            'k : [', ref%range_y%lower+1, ' : ', &
            ref%range_y%upper+1, ']'
    end if
    
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
  
  subroutine disp_box(o, prefix, ind1, ind2)
    implicit none
    type(box_info) :: o
    character(len=*), optional :: prefix
    character(len=2), optional :: ind1
    character(len=2), optional :: ind2
    character(len=2) :: ind11, ind22

    if(present(ind1)) then
       ind11 = ind1
    else
       ind11 = '2X'
    endif
    if(present(ind2)) then
       ind22 = ind2
    else
       ind22 = '6X'
    endif
    
    if(present(prefix)) then
       write(*, "("//ind11//", A)") prefix
    else
       write(*, "("//ind11//", A)") "ans =  "
    endif
    
    write(*, "("//ind22//", 3I4.1)") o%rx%lower, &
         o%ry%lower,o%rz%lower
    write(*, "("//ind22//", 3I4.1)") o%rx%upper, &
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
    res = (a%upper <= b%upper) &
         .and. (a%lower >= b%lower)
  end function
  
  !> if box a is inside box b
  function box_in(a, b) result(res)
    implicit none
    type(box_info), intent(in) :: a, b
    logical :: res

    res = (a%rx .in. b%rx) .and. &
         (a%ry .in. b%ry) .and. &
         (a%rz .in. b%rz)
    
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

  function box_num_data(a) result(res)
    implicit none
    type(box_info), intent(in) :: a
    integer :: res

    if(.not. is_valid(a)) then
       res = 0
       return
    endif

    res = (a%rx%upper - a%rx%lower + 1) * &
         (a%ry%upper - a%ry%lower + 1) * &
         (a%rz%upper - a%rz%lower + 1)
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

  subroutine disp_dist(dist, prefix)
    implicit none
    type(dist_info), intent(in) :: dist
    character(len=*), optional :: prefix

    if(present(prefix)) then
       write(*, "(A)") prefix
    else
       write(*, "(A)") "dist = "
    endif
    
    write(*, "(A, 10g2.2, A)") " lm : [", dist%lx, "]"
    write(*, "(A, 10g2.2, A)") " ln : [", dist%ly, "]"
    write(*, "(A, 10g2.2, A)") " lk : [", dist%lz, "]"
    write(*, "(A, 10g2.2, A)") "clm : [", dist%clx, "]"
    write(*, "(A, 10g2.2, A)") "cln : [", dist%cly, "]"
    write(*, "(A, 10g2.2, A)") "clk : [", dist%clz, "]"
    
  end subroutine

  function find_pos(val, idx, s1, e1) &
       result(res)
    implicit none
    integer, intent(in) :: val, s1, e1
    integer :: idx(:)
    integer :: pos
    integer :: res
    integer :: s, e

    s = s1
    e = e1
    
    pos = (s + e) / 2

    do while(.true.)
       print*, "pos1 = ", pos
       if(val >= idx(pos) .and. val < idx(pos+1)) then
          res = pos
          return
       endif

       if(pos == s .or. pos == e) then
          res = -1
          return
       endif
       
       if(val < idx(pos)) then
          e = pos
          pos = (s + e) / 2
       else
          s = pos + 1
          pos = (s + e + 1) / 2
       endif
    end do
    
  end function
  
  subroutine box_to_indices2(idx, dist, box, orange)
    implicit none
    integer, intent(out), allocatable :: idx(:,:,:)
    type(dist_info), intent(in):: dist
    type(box_info) :: box
    integer :: bx_lower, bx_upper
    integer :: by_lower, by_upper
    integer :: bz_lower, bz_upper    
    integer :: cnt, nx, ny, nz
    type(buffer_i4) :: bid_x,  bid_y,  bid_z
    type(buffer_i4) :: bid_nx, bid_ny, bid_nz
    integer :: ierr, ix, iy, iz, dim_x, dim_y, dim_z
    logical :: start
    integer :: iiz1, iiz2, iiy1, iiy2, iix1, iix2
    integer :: iix, iiy, iiz, smx, smy, smz, offset
    integer, intent(in) :: orange(:)
    integer :: i, j, k, box_shape(3)
    
    nx = size(dist%lx)
    ny = size(dist%ly)
    nz = size(dist%lz)

    if(.not. is_valid(box)) then
       allocate(idx(0,0,0))
       return
    end if
    
#:for d in ['x', 'y', 'z']
    start = .false.
    !> ${d}$-direction    
    do i = 1, size(dist%cl${d}$) - 1

       if(box%r${d}$%lower < dist%cl${d}$(i) .and. &
            box%r${d}$%upper < dist%cl${d}$(i)) then
          exit
       end if
       
       if(box%r${d}$%lower >= dist%cl${d}$(i) .and. &
            box%r${d}$%lower < dist%cl${d}$(i+1)) then
          start = .true.
          b${d}$_lower = i
       endif

       if(box%r${d}$%upper >= dist%cl${d}$(i) .and. &
            box%r${d}$%upper < dist%cl${d}$(i+1)) then
          b${d}$_upper = i
       endif
       
       if(start) then
          call push_back(bid_${d}$, &
               max(box%r${d}$%lower, dist%cl${d}$(i)))
          call push_back(bid_${d}$, &
               min(box%r${d}$%upper, dist%cl${d}$(i+1)-1))
       endif
    enddo
#:endfor

    print*, "bx1", bx_lower, "bx2", bx_upper, &
         "by1", by_lower, "by2", by_upper, &
         "bz1", bz_lower, "bz2", bz_upper
    
    dim_x = dist%clx(size(dist%clx))
    dim_y = dist%cly(size(dist%cly))
    dim_z = dist%clz(size(dist%clz))

    box_shape = shape(box)
    
    allocate(idx(box%rx%lower:box%rx%upper, &
         box%ry%lower:box%ry%upper, &
         box%rz%lower:box%rz%upper))

    cnt = 1
    i = box%rx%lower
    j = box%ry%lower
    k = box%rz%lower
    
    do iz = bz_lower, bz_upper
       do iy = by_lower, by_upper
          do ix = bx_lower, bx_upper

             iiz1 = bid_z%data((iz - bz_lower + 1) * 2 - 1)
             iiz2 = bid_z%data((iz - bz_lower + 1) * 2)
             
             iiy1 = bid_y%data((iy - by_lower + 1) * 2 - 1)
             iiy2 = bid_y%data((iy - by_lower + 1) * 2)             
             
             iix1 = bid_x%data((ix - bx_lower + 1) * 2 - 1)
             iix2 = bid_x%data((ix - bx_lower + 1) * 2)

             ! offset = iix1 + iiy1 * dim_x &
             !      + iiz1 * dim_x * dim_y

             offset = orange(ix + (iy-1) * size(dist%lx) &
                  + (iz-1)* size(dist%lx) * size(dist%ly))
             
             ! if(get_rank() == 2) then
             !    print*, "[", iix1, ",", iiy1, ",", iiz1, "]"
             !    print*, "[", iix2, ",", iiy2, ",", iiz2, "]"
             !    print*, "offset = ", offset
             ! endif
             
             ! smz = iiz2 - iiz1 + 1
             ! smy = iiy2 - iiy1 + 1
             ! smx = iix2 - iix1 + 1
             smz = dist%lz(iz)
             smy = dist%ly(iy)
             smx = dist%lx(ix)
             
             do iiz = iiz1,iiz2
                do iiy = iiy1,iiy2
                   do iix = iix1,iix2
                      
                      idx(i, j, k) = (iix - dist%clx(ix)) + &
                           (iiy - dist%cly(iy)) * smx &
                           + (iiz - dist%clz(iz)) * smx * smy &
                           + offset
                      
                      ! print*, "[", iix, ",", iiy, ",", iiz, "]"
                      ! idx(cnt) = iix + iiy * dim_x + iiz * dim_x * dim_y
                      !cnt = cnt + 1
                      i = i + 1
                   enddo
                   j = j + 1
                enddo
                k = k + 1
             enddo
          enddo
       enddo
    enddo
    
    ! !call mpi_order_start(MPI_COMM_WORLD, ierr)
    ! if(get_rank() == 2) then
    !    print*, "----------------------------------"
    !    call disp(dist, 'dist = ')
    !    call disp(box, 'box = ')
    !    print*, "bid_x = ", get_data(bid_x)
    !    print*, "bid_y = ", get_data(bid_y)
    !    print*, "bid_z = ", get_data(bid_z)
    !    print*, "idx = ", idx
    ! endif
    ! !call mpi_order_end(MPI_COMM_WORLD, ierr)
 
  end subroutine
  
end module
