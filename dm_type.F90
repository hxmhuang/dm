#include <petsc/finclude/petscmatdef.h>
#include "mat_type.h"

module dm_type 
    implicit none
    type Matrix 
		Mat  	:: x         
		integer :: xtype=MAT_XTYPE_IMPLICIT
		logical :: isGlobal=.true.
		integer :: nrow=0
		integer :: ncol=0
		integer :: ista=0
		integer :: iend=0
    end type Matrix 

end module 


