module ot_init_mod
  use ot_expr
  
contains
  subroutine ot_init(ierr)
    implicit none
    integer, intent(out) :: ierr

    call init_petsc(ierr)
    call init_kernels(ierr)    
    call init_array(ierr)
    call init_node(ierr)
    call init_expr(ierr)
    call init_geom(ierr)
  end subroutine

  subroutine ot_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr

    call finalize_data(ierr)
  end subroutine
end module
