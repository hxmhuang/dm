
#define EXIT call ot_finalize(ierr); return

module test
  type obj_ptr
     integer, pointer, dimension(:) :: p
  end type obj_ptr
  
contains
  subroutine fun1()
    print*, "hello1!"
  end subroutine
  subroutine fun2()
    print*, "hello2!"
  end subroutine
  subroutine fun3()
    print*, "hello3!"
  end subroutine
  subroutine fun4()
    print*, "hello4!"
  end subroutine

  subroutine sin()
    print*, "hello4!"
  end subroutine

  ! subroutine invoke(f)
  !   interface 
  !      subroutine f()
  !      end subroutine
  !   end interface

  !   call f()
  ! end subroutine
end module test

program main
  use ot_mod
  use test
  use ot_petsc
  use ot_test
  
  implicit none
  
  type(array) :: A, B, C, D, E, F, G, H, BB
  type(array) :: U, V, W, X, Y, Z
  integer, parameter :: NN = 1000
  integer :: ierr
  integer :: pos(2)
  type(obj_ptr), allocatable, dimension(:) :: pp
  type(array_ptr), allocatable, dimension(:) :: tp
  
  type(node) :: NA, NB, NC, ND
  integer :: v_shape(3)
  integer :: dim1, dim2, dim3, dim4
  integer :: arr1(5)
  integer :: arr(4,4)
  integer, allocatable :: aa(:)
  type(range) :: rgn
  integer(kind=4) :: ii
  character(len=100) :: aaa

  ! call invoke(transfer(loc(fun4), UNKNOWN))
  
  call ot_init(ierr)

  call test_ones()
  call test_seqs()
  call test_rand()
  
  !call test_expr()  
  !call test_math()  
  call test_slice()
  call test_set()
  !call test_find_pos()

  call ot_finalize(ierr)
  
  ! call test_size()
  ! call insert_pointer(loc(fun1), 0)
  ! call insert_pointer(loc(fun2), 1)
  ! call insert_pointer(loc(fun3), 2)
  ! call insert_pointer(loc(fun4), 3)

  ! call invoke(3)
  ! call invoke(2)
  ! call invoke(1)  
  ! call invoke(0)

  ! integer, allocatable, target :: base_r(:)
  ! integer, pointer :: r(:)=>null()

  ! r => base_r(-1:2000)
  ! print*, size(r)
  ! print*, (loc(r) - loc(base_r)) / 4
  ! print*, ((loc(r) + size(r) * 4) - loc(base_r) - 4)/4

  !pos = find_range(r(-1:2000))
  ! print*, find_range(r3(-1:2000))

  ! EXIT
  
  !I want this:
  !A(r(1:10), r(1:20)) = 10
  !A(r(1:10), r(1:20)) = B(r(100:110), r(100,120))
  !ones(2, 3) .assign. ones(2, 3)
  !call set(sub(A, r(1:10), r(1:20)), sub(B, r(1:10), r(1:20)))
  
  ! print*, find_range(r(-1:2000))
  ! print*, find_range(r(0:11))
  ! print*, find_range(r(2:15))  

  !call test_petsc_slice_dm()
  !call test_dict()
  !call test_buffer()
  !call test_ptr()
  

  !call petsc_get_shape(v_shape, A%data)
  !print*, "v_shape = ", v_shape
  ! call test_slice()
  ! EXIT

  ! call disp_info(slice(A, r(0,5), r(0,5), 0))
  
  ! B = slice(A, r(0, 5), r(0, 5), 0)
  ! call disp(B, "B = ")

  ! !tic(timer_id)
  ! call tic(1) 
  
  ! V = U + U + U;

  ! !call toc(1)
  ! call toc(1, 8*int8(NN)*NN*4, 2*int8(NN)*NN)


end program main
