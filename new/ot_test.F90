

module ot_test
  use ot_mod
  use ot_petsc
contains
  subroutine test_expr()
    implicit none
    type(tensor) :: A, B, C, D, E, F, G
    type(node) :: NA, NB, NC, ND
    integer :: m, n, k, ierr
    
    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)
    
    !call write_graph(NB, file="NB.dot")
    ! D = A + B + C
    !return

    !CALL disp('A = ', (/1,2,3,4/), orient='row')

    !three-dimensional array
    !A = ones(2, 2, 2)
    !E = 1.0/sin(ones(2,2,2) * 3.0)
    !E = ones(2,2,2) * 3.0
    !return

    NB = exp(abs(2.0 * ones(2, 2, 2) - 3.5))**2 + 1.0 / log(ones(2,2,2) * 3.0)
    call write_graph(NB, file="NB.dot")
    
    return
    
    ! C = ones(2, 2, 2)
    ! D = ones(2, 2, 2)  
    ! E = ones(2, 2, 2)

    ! call display(A, "A = ")

    C = 5.0 / (2.0 * ones(2, 2, 2))  
    call disp(C, "C = ")

    D = 10.0 * sin(ones(2, 2, 2))
    call disp(D, "D = ")


    A = ones(2, 2, 2)  
    B = ones(2, 2, 2)  
    C = ones(2, 2, 2)
    D = ones(2, 2, 2)  

    !C = 1.0 + A + B + C * D + E
    !C = 1.0 + 2. * A + 4.0 / B * C * D
    ND = 1.0 + 2. * A + 4.0 / (B * C * D) + 4 * 5 * C
    call write_graph(ND, file="ND.dot")

    !C = A - B
    !C = 2.0 + (1.0 / B) + A
    !C = 2.0 * (A - 0.1) + B  
    !E = A * B * C + (D)

    NA = B + sin(cos(exp(D+C)) + A) * 2 - 3
    NB = 2.0 + (1.0 / B) + A + tan(NA)

    call write_graph(NB, file="NB.dot")

    F = NB
    call disp(F, "F = ")

    !call write_opt_graph(NB, file="opt_graph.dot")

  end subroutine

  subroutine test_seqs()
    implicit none
    type(tensor) :: A, B, C, D
    integer :: ierr
    integer :: rank
    integer :: m, n, k

    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)
    
    rank = get_rank()
    if(rank == 0) &
         print*, '===========Test seqs=========='

    ! A = seqs(4)
    ! B = seqs(4, 3)
    C = seqs(m, n, k)
    call disp(C, 'C = ')
    ! call disp(B, 'B = ')
    ! call disp(C, 'C = ')

    ! call destroy(A, ierr)
    ! call destroy(B, ierr)
    ! call destroy(C, ierr)
  end subroutine

  subroutine test_ones()
    implicit none
    type(tensor) :: A, B, C, D
    integer :: ierr
    integer :: rank
    rank = get_rank()
    if(rank == 0) &
         print*, '===========Test ones=========='

    A = ones(4)
    B = ones(4, 3)
    C = ones(4, 4, 2)
    call disp(A, 'A = ')
    call disp(B, 'B = ')
    call disp(C, 'C = ')

    call destroy(A, ierr)
    call destroy(B, ierr)
    call destroy(C, ierr)
  end subroutine

  subroutine test_math()
    implicit none
    type(tensor) :: A, B, C, D, E
    integer :: m, n, k
    integer :: ierr
    
    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)

    
    A = seqs(m * n)
    call disp(A * 3, "A * 3 = ")
    
    ! B = A * 2
    ! call disp(B, "B = ")
    
    !call disp(3 * A, "3 * A = ")    

    ! call write_graph(A * 3, file="A3.dot")

    ! E = exp(abs(2.0 * ones(m, n, k) - 3.5))**2 + &
    !      1.0 / log(ones(m,n,k) * 3.0) + seqs(m, n, k)
    ! call disp(E, 'E = ')
    
    ! C = seqs(m, n)
    ! D = C
    ! E = C + D
    ! call disp(E, "E = ")
    
    ! C = C - D * 2
    ! call disp(C, "C = ")
    call destroy(A, ierr)
  end subroutine
  
  subroutine test_slice()
    implicit none
    type(tensor) :: A, B, C, D, E, F
    integer :: rank
    integer :: m, n, k, ierr
    
    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)

    if(get_rank() == 0) &
         print*, '===========Test slice=========='

    A = seqs(10)
    call disp(A, 'A = ')

    B = slice(A, r(1, 5))
    call disp(B, 'B = ')

    C = seqs(10, 10)
    D = slice(C, r(1, 5), r(1, 7), 1)
    call disp(D, 'D = ')

    E = seqs(10, 10, 3)
    F = slice(E, r(1, 5), r(1, 7), r(1, 2))
    call disp(F, 'F = ')

    call destroy(A, ierr)
    call destroy(B, ierr)
    call destroy(C, ierr)
    call destroy(D, ierr)
    call destroy(E, ierr)
    call destroy(F, ierr)        
  end subroutine

  subroutine test_set()
    implicit none
    type(tensor) :: A, B, C, D, E, F
    integer :: ierr, i
    integer :: rank
    integer :: m, n, k
    
    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)

    if(get_rank() == 0) &
         print*, '===========Test set=========='

    ! A = seqs(10)
    ! call disp(A, 'A = ')
    ! call set(slice(A, r(1,3)), 3)
    ! call disp(A, 'A = ')

    ! B = seqs(10, 10)
    ! call disp(B, 'B = ')
    ! call set(slice(B, r(1,3), r(1,3)), 2.1)
    ! call disp(B, 'B = ')

    ! C = seqs(5, 5, 2)
    ! call disp(C, 'C = ')
    ! call set(slice(C, r(1,3), range_all, 1), 2.1)
    ! call disp(C, 'C = ')

    ! D = ones(5, 5, 2)
    ! call set(slice(D, r(1,3),range_all,1), slice(C,r(2,4),range_all,1))
    ! call disp(D, 'D(1:3,:,1) = C(1:3,:,1)')

    ! call set(slice(D, 1, range_all, 1), slice(D, 2, range_all, 1))
    ! call disp(D, 'D = ')

    ! call set(slice(D, range_all, range_all, 2), &
    !      slice(D, range_all, range_all, 1))
    ! call disp(D, 'D = ')

    E = seqs(4,4,10)
        
    !A = slice(E, 1, 'A','A')
    !B = slice(E, 2,'A','A')
    call set(slice(E, 1, 'A', 'A'), slice(E, 2, 'A', 'A'))
    call disp(E, 'E = ')
    !call disp(B, 'B = ')
    
    return

    do i = 1, 3
       
       ! call set(slice(E, range_all, range_all, i+1), &
       !      slice(E, range_all, range_all, i))
       
       call set(slice(E, i+1,'A','A'), slice(E, i,'A','A'))

           !call MPI_Barrier(MPI_COMM_WORLD, ierr)       
    end do
    call disp(E, 'E = ')
    
    call destroy(A, ierr)
    call destroy(B, ierr)
    call destroy(C, ierr)
  end subroutine
end module
