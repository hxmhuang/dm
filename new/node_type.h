!!!file-type:f90
!the above line will help my emacs to identify this is a fortran source file

#ifndef __DM_TYPE_H__
#define __DM_TYPE_H__

#if defined(__LP64__) || defined(_LP64)
#define C_POINTER integer(8)
#else
#define C_POINTER integer(4)  
#endif

#define DEBUG
  
!define some enum types and constants  
!// integer, parameter :: type_data  = 0
!//   integer, parameter :: type_plus  = 1
!//   integer, parameter :: type_minus = 2
!//   integer, parameter :: type_mult  = 3
!//   integer, parameter :: type_divd  = 4

!here we define the node category and subcategory
! #define abstract_data   0
! #define arthimic        1
! #define math_function   2

! #define type_data   abstract_data*256+0
! #define type_ref    abstract_data*256+1
! #define type_scalar abstract_data*256+2
  
! #define type_plus   arthimic*256+0 
! #define type_minus  arthimic*256+1
! #define type_mult   arthimic*256+2
! #define type_divd   arthimic*256+3
  
! #define type_pow    math_function*256+0
! #define type_exp    math_function*256+1
! #define type_sin    math_function*256+2
! #define type_cos    math_function*256+3
! #define type_rcp    math_function*256+4

#:set i = 0  
#:include "type_def.fypp"
#:for i in range(len(L))
#define  ${L[i][0]}$ ${'\t'}$  ${L[i][1]}$
#:endfor

#define eval_tensors uplus_tensors
  
#endif
