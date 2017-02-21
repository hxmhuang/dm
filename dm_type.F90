#include <petsc/finclude/petscmatdef.h>
#include "mat_type.h"

module dm_type 
  implicit none
  type Matrix 
     Mat			:: x         
     integer			:: xtype=MAT_XTYPE_IMPLICIT
     logical 			:: isGlobal=.true.
     integer 			:: nrow=0
     integer 			:: ncol=0
     integer 			:: ista=0
     integer 			:: iend=0
     integer 			:: nx=0
     integer 			:: ny=0
     integer 			:: nz=0
  end type Matrix

  type tensor
     type(tensor), pointer :: left, right
     Vec    :: data_v
     Mat    :: data_m
     integer :: m_dim(3)
     integer :: m_shape(3)
     real(kind=8), pointer :: data1d(:), data2d(:,:), data3d(:,:,:)
  end type tensor
  
end module dm_type


