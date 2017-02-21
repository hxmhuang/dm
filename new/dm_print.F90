module dm_print_mod
  implicit none
  
  interface dm_print
     module procedure dm_print1d_i !integer 
     module procedure dm_print1d_f !float
     module procedure dm_print1d_d !double

     module procedure dm_print2d_i !integer 
     module procedure dm_print2d_f !float
     module procedure dm_print2d_d !double
     
     module procedure dm_print3d_i !integer 
     module procedure dm_print3d_f !float
     module procedure dm_print3d_d !double    
  end interface dm_print

contains
  subroutine dm_print1d_i(msg, arr)
    implicit none
    character(len=*), intent(in) :: msg
    integer,intent(in) :: arr(:)
    
    print*, msg, arr
  end subroutine

  subroutine dm_print1d_f(msg, arr)
    implicit none
    character(len=*), intent(in) :: msg
    real,intent(in) :: arr(:)
    
    print*, msg, arr

  end subroutine

  subroutine dm_print1d_d(msg, arr)
    implicit none
    character(len=*), intent(in)  :: msg
    real(kind=8),intent(in) :: arr(:)
    
    write(*,"(A)"), msg
    write(*, "(5X, 100g10.5)"), arr

  end subroutine

  subroutine dm_print2d_i(msg, arr)
    implicit none
    character(len=*), intent(in) :: msg
    integer,intent(in) :: arr(:,:)
    
    print*, msg, arr
  end subroutine

  subroutine dm_print2d_f(msg, arr)
    implicit none
    character(len=*), intent(in) :: msg
    real,intent(in) :: arr(:,:)
    integer :: i
    integer :: dim(2)

    dim = shape(arr)
    
    print*, msg

    do i = 1, dim(2)
       write(*,"(5X, 100g10.5)") arr(i, :)
    enddo

  end subroutine

  subroutine dm_print2d_d(msg, arr)
    implicit none
    character(len=*), intent(in)  :: msg
    real(kind=8),intent(in) :: arr(:,:)
    integer :: i
    integer :: dim(2)

    dim = shape(arr)
    write(*,"(2X, A)"), msg
    do i = 1, dim(2)
       write(*,"(2X, 100g15.6)") arr(i, :)
    enddo
  end subroutine
  
  subroutine dm_print3d_i(msg, arr)
    implicit none
    character(len=*), intent(in)  :: msg
    integer,intent(in), pointer :: arr(:,:,:)

    integer :: dim(3), i, j, k

    dim = shape(arr)
    
    print*, msg
    do k = 1, dim(3)
       write(*,"(4X, A, I0.1)"), "k = ",k
       do j = 1, dim(2)
          write(*,"(5X, 100g10.5)") arr(:,j,k)!)(/(arr(i,j,k), i=0,dim(1)-1)/)
       enddo
    enddo
  end subroutine

  subroutine dm_print3d_f(msg, arr)
    implicit none
    character(len=*), intent(in)  :: msg
    real,intent(in), pointer :: arr(:,:,:)

    integer :: dim(3), i, j, k

    dim = shape(arr)
    
    print*, msg
    do k = 1, dim(3)
       write(*,"(4X, A, I0.1)"), "k = ",k
       do j = 1, dim(2)
          write(*,"(5X, 100g10.5)") (/(arr(i,j,k), i=0,dim(1)-1)/)
       enddo
    enddo
  end subroutine

  subroutine dm_print3d_d(msg, arr)
    implicit none
    character(len=*), intent(in)  :: msg
    real(kind=8),intent(in), pointer :: arr(:,:,:)
    integer :: dim(3), i, j, k

    !write(*, "(Z16.6)"), loc(arr)    
    dim = shape(arr)
    
    print*, msg
    do k = 1, dim(3)
       write(*,"(4X, A, I0.1)"), "k = ",k
       do j = 1, dim(2)
          write(*,"(5X, 100g10.5)") arr(:,j,k)
       enddo
    enddo
  end subroutine
  
end module


