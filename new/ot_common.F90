module ot_common
  use mpi
  use ot_print

  
#define DEBUG
  
  public :: assert, abort, tic, toc, init
  
  interface
     subroutine abort() bind(C, name="abort")
     end subroutine abort
  end interface

  interface
     subroutine gettimeofday(ierr, val) &
          bind(C, name='gettimeofday')
       integer :: ierr
       integer*8 :: val
     end subroutine 
  end interface
  integer*8 :: timer(1000)

  interface toc
     module procedure toc1, toc2
  end interface toc

  integer :: global_counter = 0

contains

  ! -----------------------------------------------------------------------
  ! Get the input paramenters 
  ! -----------------------------------------------------------------------
  subroutine ot_option_int(str,input,ierr)
    implicit none
#include "petsc.h"
    character(len=*),intent(in) ::  str
    integer,intent(out)			::  input 
    integer,intent(out)			::  ierr 
    input=0

    call PetscOptionsGetInt(PETSC_NULL_OBJECT, &
         PETSC_NULL_CHARACTER,str,input,&
         PETSC_NULL_BOOL,ierr)
    
  end subroutine
  
  function get_global_id() result (res)
    implicit none
    integer res 
    res = global_counter
    global_counter = global_counter + 1
  end function

  function get_rank(comm) result(res)
    implicit none
    integer :: ierr
    integer :: res
    integer, optional :: comm

    if(present(comm)) then
       call MPI_Comm_rank(comm, res, ierr)
    else
       call MPI_Comm_rank(MPI_COMM_WORLD, res, ierr)       
    endif
  end function

  function get_size(comm) result(res)
    implicit none
    integer :: ierr
    integer :: res
    integer, optional :: comm
    if(present(comm)) then
       call MPI_Comm_size(comm, res, ierr)
    else
       call MPI_Comm_size(MPI_COMM_WORLD, res, ierr)
    end if
  end function

  subroutine mpi_order_start(comm, ierr)
    implicit none
    integer :: comm
    integer ,intent(out) :: ierr
    integer :: rank
    integer :: size
    integer :: i

    size = get_size(comm)
    rank = get_rank(comm)

    do i = 1, rank
       call MPI_Barrier(comm, ierr)
    end do
  end subroutine

  subroutine mpi_order_end(comm, ierr)
    implicit none
    integer :: comm
    integer ,intent(out) :: ierr
    integer :: rank
    integer :: size
    integer :: i

    size = get_size(comm)
    rank = get_rank(comm)

    do i = rank + 1, size - 1
       call MPI_Barrier(comm, ierr)
    end do
  end subroutine
  
  
  subroutine reset_global_id()
    implicit none
    global_counter = 0
  end subroutine
  
  subroutine assert(condition, file, linenum, msg)
    implicit none
    logical :: condition
    integer, optional :: linenum
    character(len=*), optional :: msg
    character(len=*), optional :: file
    
    if(.not. condition) then
       write(*,"(A)", advance="no") "Error: assertation failed"
       if(present(linenum)) then
          write(*,"(A, A, A, 4I4)", advance="no"), &
               " in ", file, " at line ",linenum
       endif
       if(present(msg)) write(*,*), "msg:", msg
       write(*,*) "."
       call abort();
    endif
  end subroutine assert
  
  subroutine tic(i) 
    integer, optional :: i
    integer :: ierr
    integer*8 :: val
    call gettimeofday(ierr, val)
    if(present(i)) then
       timer(i) = val
    else
       timer(1) = val
    endif
  end subroutine tic

  subroutine toc1(i) 
    integer,optional :: i
    integer :: ierr
    integer*8 :: val, diff

    call gettimeofday(ierr, val)
    
    if(present(i)) then
       diff = val - timer(i)
    else
       diff = val - timer(1)
    endif
    
    write(*, *) ""
    write(*, "(A, F10.6)") "Time elapsed (s) : ", real(diff,8) / 1E6
    write(*, *) ""    
  end subroutine

  subroutine toc2(i, bytes, flops) 
    integer, intent(in) :: i
    integer*8, intent(in) :: bytes, flops
    integer :: ierr
    integer*8 :: val
    real(8) :: diff
    
    call gettimeofday(ierr, val)

    !convert to seconds
    diff = real(val - timer(i), 8) / 1E6

    write(*, *) ""
    write(*, "(A20, F10.6)") "Time elapsed (s) : ", diff
    write(*, "(A20, F10.6)") "GFLOPS : ", (real(flops,8)  / diff / 1E9)
    write(*, "(A20, F10.6)") "Bandwidth (GB/s) : ", bytes / diff / 1E9
    write(*,*) ""
  end subroutine

  !> if tensor shape is [2x3x1], the function returns 2
  function find_dim(m_shape) result(dim)
    implicit none
    integer, intent(in) :: m_shape(3)
    integer :: dim
    integer :: i

    do i = size(m_shape), 1, -1
       if(m_shape(i) /= 1) then
          dim = i
          return
       end if
    end do
    dim = 0
  end function
  
  subroutine disp_shape(s, indent)
    implicit none
    integer :: s(3)
    character(len=*), optional :: indent
    character(len=10) :: header
    integer :: dim, i
    
    if(present(indent)) then
       write(header, *) "(",indent,", A)"
    else
       write(header, *) "(",indent,", A)"
    endif
    
    dim = find_dim(s)
    write(*, header, advance="no") "shape : ["
    if(dim == 0) then
       write(*, "(I0.1)", advance="no") 0
    else
       do i = 1, dim
          write(*, "(I0.1)", advance="no") s(i)
          if(i < dim) &
               write(*, "(A)", advance="no") "x"
       enddo
    endif
    write(*, "(A)") "]"
  end subroutine
  
end module ot_common
