
module ot_vector
  use ot_type
  use ot_common
  
#:for op in ['push_back', 'pop', 'remove']  
  interface ${op}$
#:for t in ['tensor', 'node']
     module procedure ${op}$_${t}$
#:if op == 'push_back'
     module procedure ${op}$_${t}$_arr
#:endif
#:endfor
  end interface ${op}$
#:endfor

  interface assign_ptr
     module procedure assign_node_ptr
     module procedure assign_tensor_ptr
  end interface assign_ptr
  
contains

#:for t in ['tensor', 'node']
  !> append array of ${t}$ pointers into vector
  subroutine push_back_${t}$_arr(v, a)
    implicit none
    type(${t}$_ptr), allocatable,  intent(inout) :: v(:)
    type(${t}$_ptr), allocatable :: tmp(:)
    type(${t}$_ptr), intent(in), dimension(:) :: a
    integer :: pos, i

    pos = 1
    if(.not. allocated(v)) then
       allocate(v(size(a)))
    else
       allocate(tmp(size(v)))
       tmp = v
       
       deallocate(v)
       allocate(v(size(tmp)+size(a)))
       
       v(1:size(tmp)) = tmp
       pos = size(tmp) + 1
       deallocate(tmp)
    end if
    
    do i = pos, size(v)
       call assign_ptr(v(i)%ptr, a(i-pos+1)%ptr)
    enddo
    
  end subroutine
#:endfor
  
#:for t in ['tensor', 'node']
  !> append a pointer of ${t}$ to the end of v
  subroutine push_back_${t}$(v, a)
    implicit none
    type(${t}$_ptr), allocatable,  intent(inout) :: v(:)
    type(${t}$_ptr), allocatable :: tmp(:)
    type(${t}$), intent(in), pointer :: a
    integer :: pos

    pos = 1
    if(.not. allocated(v)) then
       allocate(v(1))
    else
       allocate(tmp(size(v)))
       tmp = v
       deallocate(v)
       allocate(v(size(tmp)+1))
       v(1:size(tmp)) = tmp
       pos = size(tmp) + 1
       deallocate(tmp)
    end if
    v(pos)%ptr => a
  end subroutine

  subroutine assign_${t}$_ptr(p, o)
    implicit none
    type(${t}$), pointer,intent(out) :: p
    type(${t}$), target, intent(in)  :: o
    p => o
    ! add reference counter
    p%ref_cnt = p%ref_cnt + 1
  end subroutine
    
  !> pop out the first element of ${t}$_ptr array
  subroutine pop_${t}$(v, a)
    implicit none
    type(${t}$_ptr), allocatable, intent(inout) :: v(:)
    type(${t}$_ptr), allocatable :: tmp(:)
    type(${t}$), pointer, intent(out) :: a

    call assert(allocated(v), __FILE__, __LINE__, &
         "pop_${t}$: vector mush be allocated")

    a => v(1)%ptr
    allocate(tmp(size(v)-1))
    tmp = v(2:size(v))
    deallocate(v)
    allocate(v(size(tmp)))
    v = tmp
    deallocate(tmp)
  end subroutine

  !> remove the element at a  
  subroutine remove_${t}$(v, a)
    implicit none
    type(${t}$_ptr), allocatable, intent(inout) :: v(:)
    type(${t}$_ptr), allocatable :: tmp(:)
    integer :: a

    call assert(allocated(v), __FILE__, __LINE__, &
         "remove_${t}$: vector mush be allocated")

    call assert(a <= size(v), __FILE__, __LINE__, &
         "remove_${t}$: index mush be smaller than size of v")

    allocate(tmp(size(v)-1))
    tmp(1:a-1) = v(1:a-1)
    tmp(a:size(tmp)) = v(a+1:size(v))
    deallocate(v)
    allocate(v(size(tmp)))
    v = tmp
    deallocate(tmp)
  end subroutine
  
#:endfor

  
end module
