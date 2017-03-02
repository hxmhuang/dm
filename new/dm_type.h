!!!file-type:f90
!the above line will help my emacs to identify this is a fortran source file

#ifndef __DM_TYPE_H__
#define __DM_TYPE_H__

!define some enum types and constants  
!// integer, parameter :: type_data  = 0
!//   integer, parameter :: type_plus  = 1
!//   integer, parameter :: type_minus = 2
!//   integer, parameter :: type_mult  = 3
!//   integer, parameter :: type_divd  = 4

!here we define the node category and subcategory
#define abstract_data   0
#define arthimic        1
#define math_function   2


#define type_data   abstract_data*256+0
#define type_ref    abstract_data*256+1
#define type_scalar abstract_data*256+2
  
#define type_plus   arthimic*256+0 
#define type_minus  arthimic*256+1
#define type_mult   arthimic*256+2
#define type_divd   arthimic*256+3
#define type_dot    arthimic*256+4
  
#define type_pow    math_function*256+0
#define type_exp    math_function*256+1
#define type_sin    math_function*256+2
#define type_cos    math_function*256+3
#define type_rcp    math_function*256+4
  
#endif
