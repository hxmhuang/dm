
program main
  use dm_expr
  use dm_tensor
  use dm_basics
  implicit none

  type(tensor) :: A, B, C, D, E, F, G, H
  type(tensor) :: U, V, W, X, Y, Z
  integer, parameter :: NN = 1000
  integer :: ierr
  
  call dm_init(ierr)

  !three-dimensional array
  A = ones(2, 2, 2)
  call display(A, "A=")
  
  B = ones(2, 2, 2)
  
  C = A + B
  call display(C)

  !show the result of A-B
  call display(A-B, "A-B=")
  
  !one-dimensional array
  D = ones(10)
  call display(C, "C=")

  !E = A + D !error, shape does not match
  E = ones(10)
  call display(D + E)

  !generate constants matrix 5x4
  F = consts(real(2.0,kind=8), 5, 4)
  call display(F, "F=")

  !generate constants matrix 5x4  
  G = consts(real(3.0,kind=8), 5, 4)
  call display(F * G, "F*G=")
  call display(F-F/G, "F-F/G=")
  
  !performance test
  U = ones(NN, NN)

  !tic(timer_id)
  call tic(1) 
  
  V = U + U + U;

  !call toc(1)
  call toc(1, 8*int(NN,kind=8)*NN*4, 2*int(NN,kind=8)*NN)

  call tensor_destroy(A, ierr)
  call tensor_destroy(B, ierr)
  call tensor_destroy(C, ierr)
  call tensor_destroy(D, ierr)
  call tensor_destroy(E, ierr)
  call tensor_destroy(F, ierr)
  call tensor_destroy(G, ierr)
  
  call tensor_destroy(U, ierr)
  call tensor_destroy(V, ierr)

  call dm_finalize(ierr)
end program main
