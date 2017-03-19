module ot_init_mod
  use ot_expr
  
contains
  subroutine ot_init(ierr)
    implicit none
    integer, intent(out) :: ierr

    call init_data(ierr)    
    call init_tensor(ierr)
    call init_node(ierr)
    call init_expr(ierr)

  end subroutine
end module
