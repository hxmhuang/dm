
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
end module test

program main
  use ot_mod
  use test
  use ot_petsc
  use ot_test
  
  implicit none
  
  type(tensor) :: A, B, C, D, E, F, G, H, BB
  type(tensor) :: U, V, W, X, Y, Z
  integer, parameter :: NN = 1000
  integer :: ierr
  integer :: pos(2)
  type(obj_ptr), allocatable, dimension(:) :: pp
  type(tensor_ptr), allocatable, dimension(:) :: tp
  
  type(node) :: NA, NB, NC, ND
  integer :: v_shape(3)
  integer :: dim1, dim2, dim3, dim4
  integer :: arr1(5)
  integer :: arr(4,4)
  
  call ot_init(ierr)

  arr1(1) = 1;  arr1(2) = 2;  arr1(3) = 3;  arr1(4) = 4; arr1(5) = 5;

  ! arr(1,1) = 1;  arr(2,1) = 2;  arr(3,1) = 3;  arr(4,1) = 4;
  ! arr(1,2) = 1;  arr(2,2) = 2;  arr(3,2) = 3;  arr(4,2) = 4;
  ! arr(1,3) = 1;  arr(2,3) = 2;  arr(3,3) = 3;  arr(4,3) = 4;
  ! arr(1,4) = 1;  arr(2,4) = 2;  arr(3,4) = 3;  arr(4,4) = 4;
  
  !arr1(3:5) = arr1(1:3)
  print*, "arr1=", arr1
  
  ! dim1 = find_dim((/3,1,1/))
  ! dim2 = find_dim((/3,3,1/))
  ! dim3 = find_dim((/3,3,3/))
  ! dim4 = find_dim((/1,1,1/))
  ! print*, "find_dim1 = ", dim1
  ! print*, "find_dim2 = ", dim2
  ! print*, "find_dim3 = ", dim3
  ! print*, "find_dim4 = ", dim4  
  
  ! print*, "find_dim = ", find_dim((/10,1,1/))
  ! print*, "find_dim = ", find_dim((/10,10,1/))
  ! print*, "find_dim = ", find_dim((/10,10,10/))

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

  !find_range(r(-1:2000))

  !I want this:
  !A(r(1:10), r(1:20)) = 10
  !A(r(1:10), r(1:20)) = B(r(100:110), r(100,120))
  !ones(2, 3) .assign. ones(2, 3)
  !call set(sub(A, r(1:10), r(1:20)), sub(B, r(1:10), r(1:20)))
  
  ! print*, find_range(r(-1:2000))
  ! print*, find_range(r(0:11))
  ! print*, find_range(r(2:15))  

  !call test_petsc_slice_dm()
  ! call test_ones()
  ! call test_seqs()
  ! call test_slice()
  call test_set()

  EXIT

  !C = ones(2, 2, 2)

  !call petsc_get_shape(v_shape, A%data)
  !print*, "v_shape = ", v_shape

  ! call test_slice()
  ! EXIT
  
  !A = seqs(10, 10)
  !call petsc_print(A%data)

  ! B = ones(2, 2, 2)  
  ! D = A + B + C

  ! call disp(A, "A = ")
  ! call disp(C, "C = ")

  call disp_info(slice(A, r(0,5), r(0,5), 0))
  
  B = slice(A, r(0, 5), r(0, 5), 0)
  call disp(B, "B = ")


  ! C = A + B
  ! call display(C)

  ! D = log(C)
  ! call display(D, "D=")
  
  ! !show the result of A-B
  ! call display(A-B, "A-B=")
  
  ! !one-dimensional array
  ! D = ones(10)
  ! call display(C, "C=")

  ! !E = A + D !error, shape does not match
  ! E = ones(10)
  ! call display(D + E)

  ! !generate constants matrix 5x4
  ! F = consts(real(2.0, 8), 5, 4)
  ! call display(F, "F=")

  ! !generate constants matrix 5x4  
  ! G = consts(real(3.0, 8), 5, 4)
  ! call display(F * G, "F*G=")
  ! call display(F-F/G, "F-F/G=")
  
  ! !performance test
  ! U = ones(NN, NN)

  ! !tic(timer_id)
  ! call tic(1) 
  
  ! V = U + U + U;

  ! !call toc(1)
  ! call toc(1, 8*int8(NN)*NN*4, 2*int8(NN)*NN)

  ! call tensor_destroy(A, ierr)
  ! call tensor_destroy(B, ierr)
  ! call tensor_destroy(C, ierr)
  ! call tensor_destroy(D, ierr)
  ! call tensor_destroy(E, ierr)
  ! call tensor_destroy(F, ierr)
  ! call tensor_destroy(G, ierr)
  ! call tensor_destroy(U, ierr)
  ! call tensor_destroy(V, ierr)

  call ot_finalize(ierr)
end program main
