#include "config.h"

module ot_test
  use ot_mod
  use ot_petsc
  use ot_dict
contains
  subroutine test_expr()
    implicit none
    type(array) :: A, B, C, D, E, F, G
    type(expr) :: NA, NB, NC, ND
    integer :: m, n, k, ierr
    
    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)

    if(get_rank() == 0) &
         print*, '===========Test expr=========='
    
    NB = exp(abs(2.0 * ones(2, 2, 2) - 3.5))**2 &
         + 1.0 / log(ones(2,2,2) * 3.0)
    
    call write_graph(NB, file="NB.dot")

    A = ones(2, 2, 2)  
    B = ones(2, 2, 2)  
    C = ones(2, 2, 2)
    D = ones(2, 2, 2)  

    ND = 1.0 + 2. * A + 4.0 / (B * C * D) + 4 * 5 * C
    call write_graph(ND, file="ND.dot")

    !连接2个表达式
    NA = B + sin(cos(exp(D+C)) + A) * 2 - 3
    NC = 2.0 + (1.0 / B) + A + tan(NA)

    call write_graph(NC, file="NC.dot")

    !求解表达式
    F = NB
    call disp(F, "F = ")

  end subroutine

  subroutine test_seqs()
    implicit none
    type(array) :: A, B, C, D
    integer :: ierr
    integer :: rank
    integer :: m, n, k

    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)
    
    rank = get_rank()
    if(rank == 0) &
         print*, '===========Test seqs=========='

    A = seqs(m, n, k)
    
    call disp(A, 'A = ')

    ! call destroy(A, ierr)
    ! call destroy(B, ierr)
    ! call destroy(C, ierr)
  end subroutine

  subroutine test_ones()
    implicit none
    type(array) :: A, B, C, D
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

  subroutine test_rand()
    implicit none
    type(array) :: A, B, C, D
    integer :: ierr
    integer :: rank
    rank = get_rank()
    if(rank == 0) &
         print*, '===========Test rand=========='

    A = rand(4)
    B = rand(4, 3)
    C = rand(4, 4, 2)
    call disp(A, 'A = ')
    call disp(B, 'B = ')
    call disp(C, 'C = ')

    call destroy(A, ierr)
    call destroy(B, ierr)
    call destroy(C, ierr)
  end subroutine
  
  subroutine test_math()
    implicit none
    type(array) :: A, B, C, D, E
    type(expr) :: NE
    integer :: m, n, k
    integer :: ierr
    
    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)


    if(get_rank() == 0) &
         print*, '===========Test math=========='

    
    A = seqs(m, n, k) * 10 + 1
    B = 1.0 * ones(m, n, k) + 2
    C = ones(m, n, k)

    
    call disp(A, 'A = ')
    call disp(B, 'B = ')
    call disp(C, 'C = ')
    
    ! call destroy(A, ierr)
    ! return

    D = A + B + C  + sin(A + B + C)

    call disp(D, 'D = ')

    call disp(A, 'A = ')
    ! B = A * 3
    
    ! call disp(A * 3, "A * 3 = ")
    ! call disp_info(A, 'A1 = ')
    ! call disp_info(B, 'B1 = ')
    
    ! B = A * 2
    ! call disp(B, "B = ")
    
    !call disp(3 * A, "3 * A = ")    

    ! call write_graph(A * 3, file="A3.dot")

    E = exp(abs(2.0 * ones(m, n, k) - 3.5))**2 + &
         1.0 / log(ones(m,n,k) * 3.0) + seqs(m, n, k)

    call disp(E, 'E = ')
    
    ! call write_graph(cos(sin(abs(A))) + (B * 2 + 3.0), file="NE.dot")
    
    ! E = cos(sin(abs(A)))**2 + 1.0 / (B * 2 + 3.0)
    !E = cos(sin(abs(A)))**2 + (B * 2 + 3.0)
    !E = A + (B * 2 + 3.0)    

    
    !call disp_info(rcp(A))
    !E = abs(A)
    !call disp(E, 'E = ')
    
    ! C = seqs(m, n)
    ! D = C
    ! E = C + D
    ! call disp(E, "E = ")
    
    ! C = C - D * 2
    ! call disp(C, "C = ")
    
    ! call destroy(B, ierr)

    call destroy(A, ierr)
    call destroy(B, ierr)
    call destroy(C, ierr)
    call destroy(D, ierr)
    call destroy(E, ierr)
  end subroutine
  
  subroutine test_slice()
    implicit none
    type(array) :: A, B, C, D, E, F
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
    D = slice(C, r(1, 5), r(1, 7))
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
    type(array) :: A, B, C, D, E, F
    integer :: ierr, i
    integer :: rank
    integer :: m, n, k
    
    call ot_option_int('-m', m, ierr)
    call ot_option_int('-n', n, ierr)
    call ot_option_int('-k', k, ierr)

    if(get_rank() == 0) &
         print*, '===========Test set=========='

    A = seqs(m * n * k)
    call disp(A, 'A = ')
    
    A = 3
    call disp(A, 'A = ')

    B = ones(m, n, k)
    call disp(B, 'B = ')
    
    B = 2.0
    call disp(B, 'B = ')

    call set(B, 3.0)
    call disp(B, 'B = ')

    C = seqs(m, n, k)
    
    B = C + 1.0
    call disp(B, 'B = ')
    
    call set(B, 2.0 * C + 3.0)
    call disp(B, 'B = ')
    
    call set(slice(A, r(1,3)), 4)
    call disp(A, 'set(slice(A, r(1,3)), 3) = ')
    
    B = seqs(10, 10)
    call disp(B, 'B = ')
    
    call set(slice(B, r(1,3), r(1,3)), 2.1)
    call disp(B, 'B = ')

    C = seqs(5, 5, 2)
    call set(slice(C, r(1,3), 'A', 1), 2.1)
    call disp(C, 'C = ')

    ! call disp_info(C, 'info(C) = ')
    
    ! D = ones(5, 5, 2)
    ! call set(slice(D, r(1,3), 'A',1), slice(C,r(2,4), 'A', 1))
    ! call disp(D, 'D(1:3,:,1) = C(2:4,:,1)')

    ! ! call set(slice(D, 1, 'A', 1), slice(D, 2, 'A', 1))
    ! ! call disp(D, 'D = ')

    ! ! call set(slice(D, 'A', 'A', 2),  slice(D, 'A', 'A', 1))
    ! ! call disp(D, 'D = ')

    ! E = seqs(4, 4, 4)
    ! call disp(E, 'E = ')    
    ! do i = 1, 3
    !   !call set(slice(E, 'A','A', i+1), slice(E, 'A',ROW_ALL,i))       
    !    call set(slice(E, 'A','A', i+1), slice(E, 'A','A',i))
    ! end do
    ! call disp(E, 'E = ')

    ! do i = 1, 3
    !    call set(slice(E, 'A', 'A', i+1), slice(E, 'A', 'A', i))
    ! end do
    ! call disp(E, 'E = ')
    
    call destroy(A, ierr)
    call destroy(B, ierr)
    call destroy(C, ierr)
    call destroy(D, ierr)
    call destroy(E, ierr)
  end subroutine

  subroutine fun(o)
    implicit none
    integer, pointer, intent(in) :: o(:)
    integer, pointer :: o1(:)

    o1 => o
    o1(1) = o1(1) + 1

    !if you comment out this line
    !it will cause memory leak
    deallocate(o1)
  end subroutine fun
  
  subroutine test_ptr()
    implicit none
    integer, pointer :: a(:);
    integer i

    do i = 1, 1e8
       allocate(a(1000))
       a(1) = a(1) + 1
       call fun(a)
    enddo
    
    print*, "a = ", a(1)
  end subroutine test_ptr

  subroutine test_dict()
    implicit none
    type(dict_item), pointer :: A
    integer(kind=8) :: i
    
    A => null()

    do i = 0_8, 10_8
       call dict_add(A, i, int(100+i,8))
    end do

    call disp(A)
    
    ! do i = 0_8, 20_8
    !    print*, "get(i)=", dict_get(A, i)
    ! end do
  end subroutine

  subroutine test_buffer()
    use ot_buffer
    implicit none
    type(buffer_r4), pointer :: buf, buf1, buf2, buf_ptr(:), buf_ptr2
    type(buffer_list_r4) :: buf_list

    allocate(buf, buf1, buf2)
    
    ! print*, "size(buf) = ", size(buf)
    ! print*, "data = ", get_data(buf)
    
    call push_back(buf, 1.0)
    print*, "size(buf) = ", size(buf)

    call push_back(buf, 1.1)
    call push_back(buf, 1.2)
    call push_back(buf, 1.3)
    print*, "size(buf) = ", size(buf)
    print*, "data = ", get_data(buf)

    call push_back(buf1, 2.1)
    call push_back(buf1, 2.2)

    call push_back(buf2, 3.1)
    call push_back(buf2, 3.2)
    
    call push_back(buf_list, buf)
    call push_back(buf_list, buf1)
    call push_back(buf_list, buf2)

    buf_ptr2 => get_data(buf_list, 3)
    print*, "data = ", get_data(buf_ptr2)
  end subroutine

  subroutine test_box_to_indices()
    implicit none
    
  end subroutine

  subroutine test_find_pos() 
    implicit none
    integer :: idx(1000)
    integer :: i
    integer :: pos
    integer :: N

    N = 1000
    
    do i = 1,N
       idx(i) = (i - 1) * 5
    enddo

    pos = find_pos(4, idx, 1, N - 1)
    write(*, *) "pos = ", pos
    pos = find_pos(5, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(6, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(7, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(1, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(44, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(45, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(46, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(1005, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(1000, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(4995, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
    pos = find_pos(5000, idx, 1, N - 1)
    write(*, *) "pos = ", pos    
  end subroutine
end module
