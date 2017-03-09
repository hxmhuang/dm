module dm_common
  use dm_print_mod
  
  public :: assert, abort, tic, toc, init
  
  interface
     subroutine abort() bind(C, name="abort")
     end subroutine abort
  end interface

  interface
     subroutine c_gettimeofday(ierr, val) &
          bind(C, name='c_gettimeofday')
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

  function get_global_id() result (res)
    implicit none
    integer res 
    res = global_counter
    global_counter = global_counter + 1
  end function

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
    call c_gettimeofday(ierr, val)
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

    call c_gettimeofday(ierr, val)
    
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
    
    call c_gettimeofday(ierr, val)

    !convert to seconds
    diff = real(val - timer(i), 8) / 1E6

    write(*, *) ""
    write(*, "(A20, F10.6)") "Time elapsed (s) : ", diff
    write(*, "(A20, F10.6)") "GFLOPS : ", (real(flops,8)  / diff / 1E9)
    write(*, "(A20, F10.6)") "Bandwidth (GB/s) : ", bytes / diff / 1E9
    write(*,*) ""
  end subroutine

end module dm_common
