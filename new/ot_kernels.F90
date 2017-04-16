#:set kernel_file = 'kernels.fypp'
#include "type.h"
module ot_kernels
  use ot_dict
  use ot_petsc
  use ot_type
  use ot_tensor
  use ot_buffer
  
  interface
     subroutine invoke_kernel(func, A, B, &
          ops_alpha, ops_beta, ops_args, ops_num, args_num, ierr) &
          bind(C, name = 'invoke_kernel')
       use ot_type
       use ot_buffer
#:include "petsc.h"
       C_POINTER :: func
       C_POINTER :: A
       C_POINTER :: B
       C_POINTER :: ops_alpha, ops_beta
       C_POINTER :: ops_args
       C_POINTER :: ops_num
       C_POINTER :: args_num
       C_POINTER :: ierr
     end subroutine
  end interface 
  
  type(dict_item), pointer :: dict_kernels => null()
  type(dict_item), pointer :: kernels(:)
contains

  subroutine init_kernels(ierr)
    implicit none
#include "petsc.h"
    integer,intent(out) :: ierr

#:if os.path.isfile(kernel_file)
#:set lines=io.open(kernel_file).read().replace(' ','').split('\n')
#:for i in lines
#:if i.strip() != ''
#:set ii = i.split(':')  
#:set key=ii[0]
#:set func_name = 'kernel_{0}'.format(key)
#:set func=ii[-1]
    !add function pointer ${func}$
    call dict_add(dict_kernels, ${key}$_8, loc(${func_name}$))
    !call invoke_kernel(loc(${func_name}$), &
    !     0_8, (/0_8/), (/0.1_8, 0.3_8/), (/0.2_8/), (/0.3_8/), ierr)
    !print*, "ierr = ", ierr
#:endif    
#:endfor
#:endif
    
    !show registered kernels
    if(associated(dict_kernels)) &
         call disp(dict_kernels)
  end subroutine

#:if os.path.isfile(kernel_file)
#:set lines=io.open(kernel_file).read().replace(' ','').split('\n')
#:for i in lines
#:if i.strip() != ''
#!  
#:set ii = i.split(':')  
#:set key=ii[0]
#:set func=ii[-1]
#:set cnt = func.count('A')
#:set cnt_alpha = func.count('X')
#:set cnt_beta = func.count('Y')
#:set cnt_arg = func.count('B')

  !> function ${func}$
  subroutine kernel_${key}$(A, B, ops_alpha, &
       ops_beta, ops_args, args_num, ierr)
    implicit none
#include "petsc.h"
    type(tensor), intent(inout) :: A
    type(tensor_ptr), intent(in) :: B(${cnt}$)
    PetscScalar, intent(in) :: ops_alpha(${cnt_alpha}$)
    PetscScalar, intent(in) :: ops_beta(${cnt_beta}$)
    PetscScalar, intent(in) :: ops_args(args_num)
    integer, intent(in)  :: args_num
    PetscScalar, pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze
    integer, intent(out) :: ierr
#:for c in range(1,cnt+1)
    PetscScalar, pointer :: ${'{0}{1}'.format('x', c)}$(:,:,:)
#:endfor

    ! ierr = 3
    ! print*, "alpha(1) = ", ops_alpha(1)
    ! print*, "alpha(2) = ", ops_alpha(2)    
    ! print*, "hello!!${func}$"
    ! return
    
    call get_local_arr(A%data, res)
    
#:for c in range(1,cnt+1)
    call get_local_arr(B(${c}$)%ptr%data, ${'{0}{1}'.format('x', c)}$)
#:endfor
    
#:set expr = func    
#:for c in range(1,cnt+1)
#:set expr = string.replace(expr, 'A', 'x'+str(c)+'(xs:xe,ys:ye,zs:ze)', 1)
#:endfor
#!
#:for c in range(1, cnt_alpha+1)
#:set expr = string.replace(expr, 'X', 'ops_alpha({0})'.format(c), 1)
#:endfor
#!
#:for c in range(1, cnt_beta+1)
#:set expr = string.replace(expr, 'Y', 'ops_beta({0})'.format(c), 1)
#:endfor
#!
#:for c in range(1, cnt_arg+1)
#:set expr = string.replace(expr, 'B', 'ops_args(1)'.format(c), 1)
#:endfor
#!
    res(xs:xe,ys:ye,zs:ze) = ${expr}$

#:for c in range(1,cnt+1)
    call restore_local_arr(B(${c}$)%ptr%data, ${'{0}{1}'.format('x', c)}$)
#:endfor
  end subroutine
#!  xxx${key}$xxx${expr}$xxx
#:endif
#:endfor
#:endif
end module
