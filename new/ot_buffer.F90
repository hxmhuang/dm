module ot_buffer

#:set T = [['i4', 'integer'],['i8', 'integer(kind=8)'], &
  ['r4', 'real'],['r8','real(8)']]

  interface push_back
#:for type in T     
#:for stype in ['single', 'arr']
     module procedure push_back_${type[0]}$_${stype}$
#:endfor
     module procedure push_back_list_${type[0]}$     
#:endfor     
  end interface push_back

#:for type in T
  type buffer_${type[0]}$
     integer :: size = 0
     type(${type[1]}$), dimension(:), allocatable :: data
  end type buffer_${type[0]}$

  type buffer_list_${type[0]}$
     integer :: size = 0
     type(buffer_${type[0]}$), dimension(:),allocatable :: data
  end type buffer_list_${type[0]}$
  
#:endfor
  
  interface size
#:for type in T
     module procedure size_buf_${type[0]}$
     module procedure size_buf_list_${type[0]}$     
#:endfor     
  end interface size

  interface mem_size
#:for type in T
     module procedure mem_size_buf_${type[0]}$
     module procedure mem_size_buf_list_${type[0]}$     
#:endfor     
  end interface mem_size

  interface buf_data
#:for type in T
     module procedure buf_data_${type[0]}$
     module procedure buf_data_at_${type[0]}$
     module procedure buf_list_data_${type[0]}$
     module procedure buf_list_data_at_${type[0]}$
#:endfor
  end interface buf_data
  
contains

#:for type in T
#:for stype in ['single', 'arr']
  subroutine push_back_${type[0]}$_${stype}$(buf, src)
    implicit none
    type(buffer_${type[0]}$), intent(inout) :: buf
#:if stype == 'single'
    type(${type[1]}$), intent(in) :: src
#:else
    type(${type[1]}$), intent(in) :: src(:)
#:endif
    type(${type[1]}$), dimension(:), allocatable :: tmp_data
    integer :: new_len
    
    if(.not. allocated(buf%data)) then
       allocate(buf%data(1024))
       buf%size = 0
    endif
#:if stype == 'single'
    new_len = size(buf) + 1
#:else
    new_len = size(buf) + size(src)
#:endif
    
    if(size(buf%data) < new_len) then
       allocate(tmp_data(size(buf)))
       tmp_data(1:size(buf)) = buf%data(1:size(buf))
       deallocate(buf%data)
       allocate(buf%data(max(1024+size(buf), new_len)))
       buf%data(1:size(tmp_data)) = tmp_data(1:size(tmp_data))
       deallocate(tmp_data)
    endif

#:if stype == 'single'
    buf%data(buf%size+1:new_len) = src
#:else
    buf%data(buf%size+1:new_len) = src(:)
#:endif
    buf%size = new_len
  end subroutine
  
#:endfor
#:endfor

#:for type in T
  subroutine push_back_list_${type[0]}$ (buf_list, src)
    implicit none
    type(buffer_list_${type[0]}$), intent(inout) :: buf_list
    type(buffer_${type[0]}$), pointer, intent(in) :: src
    type(buffer_${type[0]}$), allocatable :: tmp_data(:)    
    integer :: new_size, old_size
    
    if(.not. allocated(buf_list%data)) then
       allocate(buf_list%data(1024))
       buf_list%size = 0
    end if

    old_size = size(buf_list)
    new_size = size(buf_list) + 1
    
    if(mem_size(buf_list) < new_size) then
       allocate(tmp_data(old_size))
       tmp_data(1:old_size) = buf_list%data(1:old_size)
       deallocate(buf_list%data)
       allocate(buf_list%data(1024+old_size))
       buf_list%data(1:old_size) = tmp_data(1:old_size)
       deallocate(tmp_data)
    endif

    buf_list%data(new_size) = src
    buf_list%size = new_size
  end subroutine
#:endfor
  
#:for type in T
  function size_buf_${type[0]}$(buf) result(res)
    implicit none
    type(buffer_${type[0]}$),intent(in) :: buf
    integer :: res
    
    res = buf%size
  end function

  function mem_size_buf_${type[0]}$(buf) result(res)
    implicit none
    type(buffer_${type[0]}$),intent(in) :: buf
    integer :: res
    
    res = size(buf%data)
  end function
  
  function size_buf_list_${type[0]}$(buf) result(res)
    implicit none
    type(buffer_list_${type[0]}$),intent(in) :: buf
    integer :: res
    
    res = buf%size
  end function

  function mem_size_buf_list_${type[0]}$(buf) result(res)
    implicit none
    type(buffer_list_${type[0]}$),intent(in) :: buf
    integer :: res
    
    res = size(buf%data)
  end function

  function buf_data_${type[0]}$(buf) result(res)
    implicit none
    type(buffer_${type[0]}$), intent(in), target :: buf
    type(${type[1]}$), pointer :: res(:)

    if(allocated(buf%data)) then
       res => buf%data(1:size(buf))
    else
       res => null()
    endif
  end function

  function buf_data_at_${type[0]}$(buf, pos) result(res)
    implicit none
    type(buffer_${type[0]}$), intent(in), target :: buf
    integer, intent(in) :: pos
    type(${type[1]}$), pointer :: res

    if(allocated(buf%data)) then
       res => buf%data(pos)
    else
       res => null()
    endif
  end function

  function buf_list_data_${type[0]}$(buf) result(res)
    implicit none
    type(buffer_list_${type[0]}$), intent(in), target :: buf
    type(buffer_${type[0]}$), pointer :: res(:)

    if(allocated(buf%data)) then
       res => buf%data(1:size(buf))
    else
       res => null()
    endif
  end function

  function buf_list_data_at_${type[0]}$(buf, pos) result(res)
    implicit none
    type(buffer_list_${type[0]}$), intent(in), target :: buf
    integer, intent(in) :: pos
    type(buffer_${type[0]}$), pointer :: res

    if(allocated(buf%data)) then
       res => buf%data(pos)
    else
       res => null()
    endif
  end function
  
#:endfor
end module
  
