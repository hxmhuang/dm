
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
  use ot_expr
  use ot_tensor
  use mod_all_func
  use dispmodule
  use test
  use ot_common
  use ot_init_mod
  use ot_type
  implicit none
  
  type(tensor) :: A, B, C, D, E, F, G, H, BB
  type(tensor) :: U, V, W, X, Y, Z
  integer, parameter :: NN = 1000
  integer :: ierr
  integer :: pos(2)
  type(obj_ptr), allocatable, dimension(:) :: pp
  type(tensor_ptr), allocatable, dimension(:) :: tp
  
  type(node) :: NA, NB, NC, ND
  
  call ot_init(ierr)

  call test_size()
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
  

  !CALL disp('A = ', (/1,2,3,4/), orient='row')
  
  !three-dimensional array
  !A = ones(2, 2, 2)
  B = exp(abs(2.0 * ones(2, 2, 2) - 3.5))**2 + log(ones(2,2,2) * 3.0)
  ! C = ones(2, 2, 2)
  ! D = ones(2, 2, 2)  
  ! E = ones(2, 2, 2)

  ! call display(A, "A = ")

  call display(B, "B = ")  
  !C = 1.0 + A + B + C * D + E
  !C = 1.0 + 2. * A + 4.0 / B * C * D
  !D = 1.0 + 2. * A + 4.0 / (B * C * D) + 4 * 5 * C
  
  !C = A - B
  !C = 2.0 + (1.0 / B) + A
  !C = 2.0 * (A - 0.1) + B  
  !E = A * B * C + (D)

  ! NA = B + sin(cos(exp(D+C)) + A) * 2 - 3
  ! NB = 2.0 + (1.0 / B) + A + tan(NA)
  ! F = NB

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

  !call dm_finalize(ierr)
end program main
