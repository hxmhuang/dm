
module ot_geom
  type corner
     integer :: xs = 0
     integer :: xe = 0
     integer :: ys = 0
     integer :: ye = 0
     integer :: zs = 0
     integer :: ze = 0
     real(8), pointer :: ptr => null()
  end type corner

end module
