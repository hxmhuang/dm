#include <petsc/finclude/petscmatdef.h>

module dm_type 

    implicit none
    type Matrix 
        Mat  x         
    end type Matrix 

    type MatrixIm 
        Mat  x         
    end type MatrixIm 

end module 


