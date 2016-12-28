program main 
  use dm
  use dm_op
  use dm_test
  implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
  type(Matrix)    		:: A,B,C,D,E,F,G,H,II,KK 
  type(Matrix)    		:: X,Y,Z,U,V,W

  integer         		:: myrank, mysize 
  integer         		:: m,n,k
  real(kind=8)    		:: ep,alpha
  real(kind=8)    		:: a1,a2,a3
  logical         		:: debug 
  integer         		:: ierr
  real(kind=8),allocatable      :: array(:)
  debug=.false.

  call dm_init1(ierr)

  call dm_comm_rank(myrank,ierr)
  call dm_comm_size(mysize,ierr)

  call dm_option_int('-m',m,ierr)
  call dm_option_int('-n',n,ierr)
  call dm_option_int('-k',k,ierr)
  call dm_option_real('-ep',ep,ierr)
  call dm_option_bool('-debug',debug,ierr)
  
  call dm_init2(m,n,k,ierr)
  if(myrank==0) then 
     print *, "==============Input paramenters==========="
     print *, "m=",m,",n=",n,",k=",k,",ep=",ep,",debug=",debug
  endif
  
  call test_dm_zeros()
  call test_dm_ones()
  call test_dm_eye()
  call test_dm_copy()
  call test_dm_add()
  call test_dm_minus()
  call test_dm_xjoin()
  call test_dm_yjoin()
  call test_dm_zjoin()
  call test_dm_mult()
  call test_dm_emult()
  call test_dm_ediv()
  call test_dm_rep()
  call test_dm_seqs()
  call test_dm_rand()
  call test_dm_sum()
  call test_dm_axpy()
  call test_dm_aypx()
  call test_dm_trans()
  call test_dm_xyt()
  call test_dm_xty()
  call test_dm_exp()
  call test_dm_log()
  call test_dm_sqrt()
  call test_dm_squ()
  call test_dm_cube()
  call test_dm_solve()
  call test_dm_setvalue()
  call test_dm_getsub()
  call test_dm_setvalues()
  call test_dm_getvalues()
  call test_dm_norm()
  call test_dm_lt()
  call test_dm_le()
  call test_dm_gt()
  call test_dm_ge()
  call test_dm_eq()
  call test_dm_nq()
  call test_dm_sparse()
  ! call test_dm_save()
  ! call test_dm_load()
  ! call test_dm_save3d()
  ! call test_dm_load3d()

  call test_OP_AXF()
  call test_OP_AXB()
  
  call test_OP_AYF()
  call test_OP_AYB()
  
  !call test_OP_AZF()
  !call test_OP_AZB()

  call test_OP_DXF()  
  call test_OP_DXB()

  call test_OP_DYF()  
  call test_OP_DYB()
  
  ! call test_OP_DZB()  
  ! call test_OP_DZF()

  print *,""

  if(myrank==0) then
     print *, "****************************************"     
     print *, "*                                      *"
     print *, "*            Test Operator Module      *"
     print *, "*                                      *"
     print *, "****************************************"
  endif
  
  call CreateOperators(2*m+1, 2*n+1, k)

  call test_AXF()
  call test_AXB()

  call test_AYF()
  call test_AYB()

  ! call test_AZF()
  ! call test_AZB()
  
  call test_DXF()
  call test_DXB()
  call test_DXC()
  
  call test_DYF()
  call test_DYB()
  call test_DYC()
  
  ! call test_DZF()
  ! call test_DZB()
  ! call test_DZC()

  
  ! if(myrank==0) print *, "==============Test dm_cart2sph============"
  ! filename="md001.00004"
  ! call dm_load(filename,.true., A, ierr)	
  ! call dm_load(filename,.false.,C, ierr)
  
  ! call dm_cart2sph(A,B,ierr)
  ! call dm_cart2sph(C,D,ierr)
  
  ! if(debug) then
  !    if(myrank==0) print *, ">A=dm_load(filename,A,.true.ierr) "
  !    call dm_view(A,ierr)
  !    if(myrank==0) print *, ">B=dm_cart2sph(A)"
  !    call dm_view(B,ierr)
  !    if(myrank==0) print *, ">C=dm_load(filename,A,.false.ierr) "
  !    if(myrank==0) call dm_view(C,ierr)
  !    if(myrank==0) print *, ">D=dm_cart2sph(A)"
  !    if(myrank==0) call dm_view(D,ierr)
  ! endif
  
  ! call dm_destroy(A,ierr)
  ! call dm_destroy(B,ierr)
  ! call dm_destroy(C,ierr)
  ! call dm_destroy(D,ierr)

  
  call dm_finalize(ierr)
end program main
