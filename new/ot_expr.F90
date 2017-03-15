module ot_expr
  use ot_common
  use ot_tensor
  use dispmodule
  use ot_type
  use ot_vector
  
#include "node_type.h"
!#define ASSERT(x,msg) call assert(x, __FILE__, __LINE__, msg)
#:set DEBUG = 1
#define DEBUG
#:include "type_def.fypp"
  
  public :: operator(+),operator(-)
  public :: operator(*),operator(/)
  public :: assignment(=)

  interface display
     module procedure display1, display2
  end interface
  
  interface assignment (=)
     module procedure tensor_assign_tensor
     module procedure node_assign_tensor
     module procedure node_assign_node
  end interface assignment (=)

  interface node_new
     module procedure node_new1, node_new2, node_new3
     module procedure node_new4, node_new5, node_new6
     module procedure node_new7, node_new8     
  end interface node_new

  interface reallocate
     module procedure reallocate1,reallocate2
  end interface reallocate

  interface disp_info
     module procedure disp_info_tensor
     module procedure disp_info_node     
  end interface disp_info

  ! interface operator (-)
  !    module procedure node_uminus, tensor_uminus
  ! end interface operator (-)

  ! interface operator (+)
  !    module procedure node_uplus, tensor_uplus
  ! end interface operator (+)

  interface operator (**)
#:for type1 in ['node', 'tensor']     
#:for type2 in ['integer', 'real', 'real8']
     module procedure ${type1}$_pow_${type2}$
#:endfor
#:endfor
  end interface operator (**)
  
  !the following code using preprossor to create interfaces
#:for op in [['plus','+'], ['minus','-'], ['mult','*'], ['divd','/'], &
  ['gt','>'], ['ge', '>='], ['lt', '<'],['le', '<='],['eq','=='],['ne','/=']]
  interface operator (${op[1]}$)
#:for type1 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:for type2 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:if any(type1 in s for s in ['real(8)', 'real', 'integer']) & 
     and any(type2 in s for s in ['real(8)', 'real', 'integer'])
#:else
#:set name1 = re.sub('\(8\)', '8', type1)
#:set name2 = re.sub('\(8\)', '8', type2)
     module procedure ${name1}$_${op[0]}$_${name2}$
#:endif
#:endfor
#:endfor
  end interface operator (${op[1]}$)
  
#:endfor

!interface for exp, log, pow, ...
#:for e in L
#:if e[1] >= 111
#:if e[3] == '+' or e[3] == '-'
  interface operator(${e[3]}$)
#:else
  interface ${e[3]}$
#:endif
#:for type1 in ['tensor', 'node']
     module procedure ${e[2]}$_${type1}$
#:endfor
  end interface

#:endif
#:endfor

  character(len=8), dimension(${L[len(L)-1][1]}$) :: op_names
  
contains

  subroutine init_expr(ierr)
    use ot_tensor
    implicit none
    integer, intent(out) :: ierr

#:for e in L
    op_names(${e[1]}$) = "${e[3]}$"
#:endfor
    
    call init_tensor(ierr)
    
  end subroutine
  
  subroutine node_new1(B, A) 
    implicit none
    type(tensor), intent(in), target :: A
    type(node),   intent(out) :: B
    B%m_dim = A%m_dim
    B%m_shape = A%m_shape
    B%data => A
    B%local_block = A%local_block
    B%node_type = type_data
  end subroutine

  !copy node A to node B
  subroutine node_new2(B, A) 
    implicit none
    type(node), intent(in)   :: A
    type(node), intent(out)  :: B
    integer i
    
    B%m_dim = A%m_dim
    B%m_shape = A%m_shape
    B%alpha = A%alpha
    B%beta = A%beta
    B%data => A%data
    B%scalar = A%scalar
    B%args = A%args
    B%node_type = A%node_type
    B%local_block = A%local_block
    
    call copy_node_ptr(B%opt_operands, A%opt_operands)
    call copy_node_ptr(B%operands, A%operands)
  end subroutine

  subroutine node_new3(C, op_type, left, right)
    implicit none
    integer, intent(in) :: op_type
    type(node), intent(in), target :: left, right
    type(node), intent(out) :: C
    real(8) :: x

    if((.not. is_scalar(left)) .and. &
         (.not. is_scalar(right))) then
       if(left%m_dim /= right%m_dim .or. &
            left%m_dim /= right%m_dim) then
          
          write(*,*) "left%type=", left%node_type, &
               "right%type=", right%node_type
          write(*,*) "left%dim=", left%m_shape, &
               "right%dim=", right%m_shape
          
       end if
       
       call assert(left%m_dim == right%m_dim .and. &
            left%m_dim == right%m_dim, &
            __FILE__, __LINE__, &
            "shape of the left and right leaf does not match.")
    endif
    
    if(is_scalar(left)) then
       x = left % scalar
       select case(op_type)
       case (type_plus) ! x + A
          call node_new(C, right)
          C%alpha = x + C%alpha
       case (type_minus) ! x - A
          call node_new(C, right)
          C%alpha = x - C%alpha
          C%beta  = -C%beta
       case (type_mult) ! x .* A
          call node_new(C, right)
          C%alpha = x * C%alpha
          C%beta  = x * C%beta
       case (type_divd) ! x ./ A
          C%m_dim   = right%m_dim
          C%m_shape = right%m_shape
          C%is_implicit = .true.
          C%beta = x
          C%node_type = type_rcp
          ! allocate(C%operands(1))
          ! C%operands(1)%ptr => right
          call push_back(C%operands, right)
       end select
    else if(is_scalar(right)) then
       x = right % scalar
       call node_new(C, left)      
       select case(op_type)
       case (type_plus) ! A + x
          C%alpha = C%alpha + x
          C%beta  = C%beta
       case (type_minus) ! A - a
          C%alpha = C%alpha - x
          C%beta  = C%beta
       case (type_mult) ! A .* a
          C%alpha = C%alpha * x
          C%beta  = C%beta  * x
       case (type_divd)
          C%alpha = C%alpha / x
          C%alpha = C%beta  / x
       end select
    else
       allocate(C%operands(2))
       C%operands(1)%ptr => left
       C%operands(2)%ptr => right

       if(is_scalar(left)) then
          C%m_dim   = right%m_dim
          C%m_shape = right%m_shape
       else
          C%m_dim   = left%m_dim
          C%m_shape = left%m_shape
       endif
       C%node_type = op_type
       C%is_implicit = .true.
    end if
  end subroutine

  !> create a node from a real*8 scalar
  subroutine node_new4(B, scalar) 
    implicit none
    real(8), intent(in), target :: scalar
    type(node), intent(out) :: B
    
    B%m_dim = 0
    B%m_shape = (/1, 1, 1/)
    B%scalar = scalar
    B%node_type = type_scalar
  end subroutine

  !> create a node from a real*8 scalar  
  subroutine node_new5(B, scalar)
    implicit none
    real, intent(in):: scalar
    type(node),intent(out) :: B

    call node_new4(B, real(scalar,8))
  end subroutine

  !> create a node from a real*8 scalar  
  subroutine node_new6(B, scalar)
    implicit none
    integer, intent(in) :: scalar
    type(node), intent(out) :: B

    call node_new4(B, real(scalar,8))
  end subroutine

  subroutine node_new7(dst, op_type, src)
    implicit none
    type(node), intent(out) :: dst
    integer :: op_type
    type(tensor), intent(in) :: src
    type(node) :: dst1
    
    call node_new(dst1, src)
    
    call node_new(dst, op_type, dst1)

    !print*, trim(op_names(op_type))
  end subroutine

  subroutine node_new8(dst, op_type, src)
    implicit none
    type(node), intent(out) :: dst
    integer :: op_type
    type(node), intent(in) :: src
    type(node), pointer :: p
    
    dst%node_type = op_type
    dst%m_dim     = src%m_dim
    dst%m_shape   = src%m_shape
    dst%is_implicit = .true.
    dst%local_block = src%local_block

    allocate(p)
    call node_new(p, src)

    call push_back(dst%operands, p)    
    ! allocate(dst%operands(1))
    ! dst%operands(1)%ptr => p
  end subroutine

  subroutine disp_info_node(o, prefix)
    implicit none
    type(node), intent(in) :: o
    character(len=*), optional, intent(in) :: prefix
    integer :: i
    
    write(*,*) ""
    if(present(prefix)) then
       write(*, "(A)") prefix
    else
       write(*, "(A)") "ans = "
    endif

    write(*, "(4X, A, I0.1)") "node id : ", o%id
#ifdef DEBUG
    write(*, "(4X, A, Z16.16)") "obj addr : 0x", loc(o)
#endif
    write(*, "(4X, A, A)") "node type : ", op_names(o%node_type)    
    write(*, "(4X, A)", advance="no") "shape : ["
    do i = 1, o%m_dim
       write(*, "(I0.1)", advance="no") o%m_shape(i)
       if(i < o%m_dim) write(*, "(A)", advance="no") "x"
    enddo
    write(*, "(A)") "]"
    
    write(*, "(4X, A, L1)") "is_implicit : ", o%is_implicit

    write(*, "(4X, A, F8.4)") "alpha : ", o%alpha
    write(*, "(4X, A, F8.4)") "beta  : ", o%beta
    if(is_scalar(o)) write(*, "(A, F0.8)") "scalar  : ", o%scalar
    write(*, "(4X, A, 10F8.4)") "args : ", o%args(1:5)

    if(.not. is_data(o)) then
       write(*, "(4X, A)") "operands : "
       do i = 1, size(o%operands)
          write(*, "(6X, A, Z16.16)") "0x", loc(o%operands(1)%ptr)
       enddo
    else
       write(*, "(4X, A, Z16.16)") "data : 0x", loc(o%data)
    endif

    write(*, "(A)") ""
  end subroutine

  subroutine display1(A, msg)
    implicit none
    type(node) :: A
    character(len=*),intent(in),optional :: msg
    if(present(msg)) then
       call display2(A%data, msg)
    else
       call display2(A%data)
    end if
  end subroutine

  subroutine node_assign_node(A, B)
    implicit none
    type(node), intent(inout), target :: A
    type(node), intent(in),    target :: B

#:if DEBUG > 0
    print*, "call node_assign_node"
#:endif
    
    A%data  => B%data
    A%node_type = B%node_type
    A%m_dim = B%m_dim
    A%m_shape = B%m_shape
    A%alpha = B%alpha
    A%beta = B%beta
    A%scalar = B%scalar
    A%is_implicit = B%is_implicit
    call copy_node_ptr(A%operands, B%operands)
    call copy_node_ptr(A%opt_operands, B%opt_operands)
    
  end subroutine
  
  subroutine node_assign_tensor(A, B)
    implicit none

    type(tensor), intent(out) :: A
    
    !the type of B must be "intent(in)"
    type(node), intent(in)   :: B 
    type(node) :: C
    type(node_ptr), allocatable :: res(:)
    integer ierr
    
    call node_new(C, B)

    call node_combine(C, res)
    
    call eval(A, C, ierr)

    !call display(C%data)
    
    ! call node_destroy(C, ierr)
    
    !call tensor_copy(A, B%data)
    
    !if(B%node_type /= type_data) call eval(A)
    ! if(B%node_type == type_data) then
    !    call tensor_copy(A, B)
    ! else
    !    call tensor_copy(C, B)
    !    call eval(C)
    !    call tensor_copy(A, C)
    ! endif
    
    if(allocated(res)) deallocate(res)    
  end subroutine

  recursive subroutine eval(res, A, ierr, is_root_arg)
    implicit none
    type(tensor), intent(inout) :: res
    type(node), intent(in) :: A
    integer :: i, num_oprands
    integer, intent(out) :: ierr
    real(8), allocatable :: ops_alpha(:), ops_beta(:)
    real(8) :: alpha, beta
    type(tensor_ptr), allocatable :: ops_tensor(:)
    logical, optional, intent(in) :: is_root_arg
    logical :: is_root
    type(node), pointer :: node_ptr
    
    ! interface
    !    subroutine invoke(p, res, operands, alpha, beta, args, n) &
    !         bind(C, name="invoke")
    !      implicit none
    !      integer :: p
    !      C_POINTER :: res
    !      C_POINTER :: operands
    !      C_POINTER :: alpha
    !      C_POINTER :: beta
    !      C_POINTER :: args
    !      integer :: n
    !    end subroutine
    ! end interface

    if(present(is_root_arg)) then
       is_root = is_root_arg
    else
       is_root = .true.
    endif

    call disp_info(A, 'A = ')
    
    !if this node is the root node, force to evaluate
    if(is_data(A) .and. .not. is_root) return

    !if(allocated(A%operands)) then
    if(.not. is_data(A)) then       
       num_oprands = size(A%operands)

       allocate(ops_tensor(num_oprands), &
            ops_alpha(num_oprands), ops_beta(num_oprands))
       
       !recusively evaluate sub-nodes
       do i = 1, num_oprands

          ! call disp_info(A, 'A = ')
          ! call disp_info(A%operands(i)%ptr, 'A1 = ')
          
          node_ptr => A%operands(i)%ptr

          if(.not. associated(node_ptr%data)) then
             allocate(node_ptr%data)
          endif
          
          call eval(node_ptr%data, node_ptr, ierr, .false.)
          
          ops_tensor(i)%ptr => node_ptr%data
          ops_alpha(i) = node_ptr%alpha
          ops_beta(i)  = node_ptr%beta
       end do
    else
       allocate(ops_tensor(1), ops_alpha(1), ops_beta(1))
       ops_tensor(1)%ptr => A%data
       ops_alpha(1) = A%alpha
       ops_beta(1)  = A%beta
    endif

    call tensor_duplicate(res, ops_tensor(1)%ptr)
    
    !call invoke(A%node_type, loc(A%data), loc(ops_tensor), &
    !     loc(alpha), loc(beta), loc(A%args), num_oprands)
    !write(*, "(A, Z16.6)") "A2=", loc(A%data)

    if(is_root) then
       alpha = A%alpha
       beta = A%beta
    else
       alpha = 0
       beta = 1.0
    endif
    
    if(need_eval(A) .or. is_root) then
       select case(A%node_type)
#:for e in L
#:if e[4]=='A' or e[4]=='B' or e[4]=='C'
       case (${e[1]}$)
          call ${e[2]}$_tensors(res, &
               ops_tensor, ops_alpha, ops_beta, A%args, alpha, beta)
#:endif
#:endfor
       case default
          call eval_tensors(res, &
               ops_tensor, ops_alpha, ops_beta, A%args, alpha, beta)
       end select

       ! if(.not. associated(A%data)) then
       !    print*, "execute function ", op_names(A%node_type), " failed."
       !    call abort()
       ! end if
    endif
      
    print*, "evalation of ", trim(op_names(A%node_type)), " is complete."
  end subroutine
  
  !> this function will be called in the following two situations:
  !> 1) A = B, where A and B are both tensors
  !> 2) func(A), where A is a tensor and passed to a function as an argument.
  !> in both case, we do not need a deep copy
  subroutine tensor_assign_tensor(A, B)
    implicit none
    type(tensor), intent(out) :: A
    type(tensor), intent(in)  :: B

    print*, "--> called tensor_assign_tensor"
    
    if(B%is_implicit) then
       call tensor_copy(A, B)
    else
       call tensor_deep_copy(A, B)
    end if
  end subroutine
  
  recursive subroutine node_combine(A, res) 
    implicit none
    type(node), intent(inout), target :: A
    type(node_ptr),  allocatable, intent(out) :: res(:)
    type(node_ptr),  allocatable :: left(:), right(:)
    type(node_ptr),  allocatable :: subs(:), tmp(:)
    integer i, cnt

    cnt = 0
    
    if(A%node_type == type_data &
         .or. A%node_type==type_scalar) return

    do i = 1,size(A%operands)
       !combine the left node
       call node_combine(A%operands(i)%ptr, subs)

       ! if the node_type does not match its leaf node_type,
       ! not accept the returned list from leaf
       if(A%node_type /= A%operands(i)%ptr%node_type .or. &
            (.not. is_arithmetic(A))) then
          if(allocated(subs)) deallocate(subs)
          allocate(subs(1))
          subs(1)%ptr => A%operands(i)%ptr
       end if
       
       if(allocated(subs)) then
          cnt = cnt + size(subs)

          call reallocate(res, cnt)
          res(cnt-size(subs)+1 : cnt) = subs
       endif
    enddo

    call copy_node_ptr(A%opt_operands, res)

! #ifdef DEBUG
!     write(*,*) "**********************"
!     write(*, "(1X, Z16.16, A, I0.3)"), loc(A), " type=", A%node_type
    
!     if(allocated(A%opt_operands)) then
!        do i = 1, size(A%opt_operands)
!           write(*, "(5X, Z16.16, A, I0.3)"), &
!                loc(A%opt_operands(i)%ptr), &
!                " type=", A%opt_operands(i)%ptr%node_type 

!        enddo
!     endif
! #endif
    
  end subroutine


!   recursive subroutine node_optimize(A, res) 
!     implicit none
!     type(node), intent(inout), target :: A
!     type(node_ptr),  allocatable, intent(out) :: res(:)
!     type(node_ptr),  allocatable :: subs(:), tmp(:)
!     integer i, cnt, num_oprands
!     type(node), pointer :: child
    
!     cnt = 0
    
!     if(is_data(A)) return

!     num_oprands = size(A%operands)

!     do i = 1, num_oprands

!        child => A%operands(i)%ptr
       
!        if(A%node_type == child%node_type .and. &
!             is_arithmetic(A)) then
!           call reallocate(A%operands, num_oprands+size(A%operands))
!           A%operands()
!        end if
       
!        ! if the node_type does not match its leaf node_type,
!        ! not accept the returned list from leaf
!        if(A%node_type == A%operands(i)%ptr%node_type) then
!           if(allocated(subs)) deallocate(subs)
!           allocate(subs(1))
!           subs(1)%ptr => A%operands(i)%ptr
!        end if
       
!        if(allocated(subs)) then
!           cnt = cnt + size(subs)

!           call reallocate(res, cnt)
!           res(cnt-size(subs)+1 : cnt) = subs
!        endif
!     enddo

!     do i = 1, size(A%operands) 
!        !combine the left node
!        call node_optimize(A%operands(i)%ptr, subs)
!     enddo
!     call copy_node_ptr(A%opt_operands, res)

! #ifdef DEBUG
!     write(*,*) "**********************"
!     write(*, "(1X, Z16.16, A, I0.3)"), loc(A), " type=", A%node_type
    
!     if(allocated(A%opt_operands)) then
!        do i = 1, size(A%opt_operands)
!           write(*, "(5X, Z16.16, A, I0.3)"), &
!                loc(A%opt_operands(i)%ptr), &
!                " type=", A%opt_operands(i)%ptr%node_type 

!        enddo
!     endif
! #endif
!    end subroutine
  
  recursive subroutine node_destroy(A, ierr)
    implicit none
    type(node), intent(inout), target :: A
    integer :: ierr
    integer i

    if(A%node_type == type_data .or. &
         A%node_type == type_scalar) return
    
    ! do i = 1, size(A%opt_operands)
    !    write(*, "(Z16.16, A, I4, A, I4)"), &
    !         loc(A%opt_operands(i)%ptr), " = ", A%opt_operands(i)%ptr%node_type, &
    !         " = ", size(A%opt_operands(i)%ptr%opt_operands)
    ! enddo
    
    if(allocated(A%operands)) then
       do i = 1, size(A%operands)
          call node_destroy(A%operands(i)%ptr, ierr)
       enddo
    endif
    
    !destroy the operands
    deallocate(A%opt_operands)
    !A%opt_operands => null()

    !destroy data if the data is implicit
    if(associated(A%data)) then
       if(A%data%is_implicit) then
          call tensor_destroy(A%data, ierr)
          A%data => null()
       end if
    end if
  end subroutine

#:for o in ['', 'opt_']
  recursive subroutine write_${o}$graph(A, is_root, file)
    implicit none
    type(node), intent(inout) :: A
    integer :: i
    character(len=100) :: label1, label2
    character(len=*),optional :: file
    integer id
    logical,intent(in), optional :: is_root
    logical :: am_i_root
    integer, parameter :: out_unit = 20
    character(len=10) :: op_name
    
    character(len=*), parameter :: &
         format_label = "(I0.1,A,A,A,I0.1,A,F4.1,A,F4.1,A)"
    character(len=*), parameter :: &
         format_label1 = "(I0.1,A,A,A)"
    character(len=*), parameter :: &
         format_map = "(I0.1,A,I0.1,A)"
    
    if(present(is_root)) then
       am_i_root = is_root
    else
       am_i_root = .true.
    endif

    if(am_i_root) then
       call reset_global_id()
       if(present(file)) then
          open(unit=out_unit,file=file, &
               action="write", status="replace")
       else
          open(unit=out_unit,file="graph.dot", &
               action="write", status="replace")
       endif
       write(out_unit, *) "digraph G {"
    end if

    write(op_name, *) trim(op_names(A%node_type))

    id = get_global_id()
    A%id = id
    
#:set SHOW_ALPHA_BETA = 1
    
#:if SHOW_ALPHA_BETA > 0
    write(out_unit, format_label) id, &
         '[label="[',trim(adjustl(op_name)),']\n id=', A%id, '\n(',&
         A%alpha, ',',A%beta, ')"];'
#:else
    write(out_unit, format_label1) id, &
         '[label="',trim((adjustl(op_name)), '"];'
#:endif
    
    !print*, label1
    if(is_data(A)) then
       if(am_i_root) then
          write(out_unit, *) "}"
          close (out_unit)
       end if
       return
    end if
    
! #ifdef DEBUG
! #:if o == "opt_"
!     print*, "-------------------------"
! #:else
!     print*, "~~~~~~~~~~~~~~~~~~~~~~~~~"
! #:endif
!     do i = 1, size(A%${o}$operands)
!        write(*, "(I0.1, A, Z16.16)"), &
!             i," = ", loc(A%${o}$operands(i)%ptr)           
!     end do
! #endif

    if(allocated(A%${o}$operands)) then
       do i = 1, size(A%${o}$operands)
          call write_${o}$graph(A%${o}$operands(i)%ptr,  .false.)
          write(out_unit, format_map), &
               A%id, "->", A%${o}$operands(i)%ptr%id, ";"
       end do
    endif
    
    if(am_i_root) then
       write(out_unit, *) "}"
       close (out_unit)
    end if
    
  end subroutine
#:endfor
  
  function is_scalar(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    res = (A%node_type == type_scalar)
  end function

  function is_data(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    res = (A%node_type == type_scalar) .or. &
         (A%node_type == type_data) .or. &
         (A%node_type == type_ref) 
  end function

  function is_arithmetic(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    res = (A%node_type == type_plus) .or. &
         (A%node_type == type_minus) .or. &
         (A%node_type == type_mult)  .or. &
         (A%node_type == type_divd)
  end function

  function need_eval(A) result(res)
    implicit none
    type(node) :: A
    logical :: res
    integer :: i

    do i = 1, size(A%args)
       if(A%args(i) .ne. 0.) then
          res = .true.
          return
       endif
    enddo
    
    res = (A%alpha /= 0 .or. &
         A%beta /= 1.0 .or. &
         (.not. associated(A%data)))
    
  end function
  
  subroutine copy_node_ptr(dst, src)
    implicit none
    type(node_ptr), allocatable, intent(in) :: src(:)
    type(node_ptr), allocatable, intent(inout) :: dst(:)

    if(allocated(src)) then
       if(allocated(dst)) deallocate(dst)
       allocate(dst(size(src)))
       dst = src
    end if
  end subroutine

  !> reallocate an array while keeping its data
  subroutine reallocate1(arr, len)
    implicit none
    type(node_ptr), allocatable, intent(inout) :: arr(:)
    type(node_ptr), allocatable :: tmp(:)
    integer :: len, copy_len
    
    if(allocated(arr)) then
       !save the data to tmp
       allocate(tmp(size(arr)))
       tmp = arr
       deallocate(arr)
       
       !allocate new size
       allocate(arr(len))

       !copy data back from tmp
       copy_len = min(size(tmp), len)
       arr(1:copy_len) = tmp(1:copy_len)
       deallocate(tmp)
    else
       allocate(arr(len))
    endif
  end subroutine

  !> reallocate an array while keeping its data
  subroutine reallocate2(arr, len)
    implicit none
    type(node_ptr), pointer, intent(inout) :: arr(:)
    type(node_ptr), pointer :: tmp(:)
    integer :: len, copy_len
    
    if(associated(arr)) then
       !save the data to tmp
       allocate(tmp(size(arr)))
       tmp = arr
       deallocate(arr)
       
       !allocate new size
       allocate(arr(len))

       !copy data back from tmp
       copy_len = min(size(tmp), len)
       arr(1:copy_len) = tmp(1:copy_len)
       deallocate(tmp)
    else
       allocate(arr(len))
    endif
  end subroutine

  !the following code using preprossor to create subroutines  
#:for type1 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:for type2 in ['real(8)', 'real', 'integer', 'tensor', 'node']
#:for op in ['plus', 'minus', 'mult', 'divd', 'gt', 'ge', 'lt', 'le', 'eq', 'ne']
#:if any(type1 in s for s in ['real(8)', 'real', 'integer']) & 
     and any(type2 in s for s in ['real(8)', 'real', 'integer'])  
#:else  
#:set name1 = re.sub('\(8\)', '8', type1)
#:set name2 = re.sub('\(8\)', '8', type2)
  function ${name1}$_${op}$_${name2}$ (A, B) result(res)
    implicit none       
    type(${type1}$),intent(in),target  :: A
    type(${type2}$),intent(in),target  :: B
    type(node), target :: res
    type(node), pointer :: C, D

    allocate(C, D)

    call node_new(C, A)
    call node_new(D, B)
    call node_new(res, type_${op}$, C, D)
  end function
  
#:endif
#:endfor
#:endfor
#:endfor

#:for e in L
#:if e[1] >= 111
  #:for type in ['tensor', 'node']
  !> function ${e[2]}$(${type}$), return a tree node 
  function ${e[2]}$_${type}$ (o) result(res)
    implicit none
    type(${type}$), intent(in) :: o
    type(node) :: res

    call node_new(res, ${e[0]}$, o)
    
  end function

#:endfor
#:endif
#:endfor

#:for type1 in ['node', 'tensor']
#:for type2 in ['integer', 'real', 'real(8)']
#:set name2 = re.sub('\(8\)', '8', type2)
  !> power function ${type1}$**${type2}$
  function ${type1}$_pow_${name2}$(o, n) result(res)
    implicit none
    type(${type1}$), intent(in) :: o
    type(${type2}$), intent(in) :: n
    type(node) :: res
    type(node), pointer :: new_o

    allocate(new_o)
    call node_new(new_o, o)
    call node_new(res, type_pow, new_o)
    res%args(1) = n
    
  end function
  
#:endfor
#:endfor
end module

