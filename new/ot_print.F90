module ot_print
  implicit none
  
  interface printx
#:for t in [['int', 'integer'], ['real', 'real'], ['real8', 'real(8)']]
     module procedure ot_print1d_${t[0]}$
     module procedure ot_print2d_${t[0]}$
     module procedure ot_print3d_${t[0]}$     
#:endfor
  end interface printx

contains

#:for t in [['int', 'integer'], ['real', 'real'], ['real8', 'real(8)']]  
  subroutine ot_print1d_${t[0]}$(msg, arr)
    implicit none
    character(len=*), intent(in)  :: msg
    type(${t[1]}$),intent(in),pointer :: arr(:)
    
    write(*,"(A)") msg
    write(*, "(5X, 100g12.5)") arr

  end subroutine

  subroutine ot_print2d_${t[0]}$(msg, arr)
    implicit none
    character(len=*), intent(in)  :: msg
    type(${t[1]}$),intent(in),pointer :: arr(:,:)
    integer :: i
    integer :: dim(2)

    dim = shape(arr)
    write(*,"(2X, A)") msg
    do i = 1, dim(2)
       write(*,"(2X, 100g15.6)") arr(i, :)
    enddo
  end subroutine

  subroutine ot_print3d_${t[0]}$(msg, arr)
    implicit none
    character(len=*), intent(in)  :: msg
    type(${t[1]}$),intent(in), pointer :: arr(:,:,:)

    integer :: dim(3), i, j, k

    dim = shape(arr)

    print*, "dim=", dim
    
    print*, msg
    do k = 1, dim(3)
       write(*,"(4X, A, I0.1)") "k = ",k
       do j = 1, dim(2)
          write(*,"(5X, 100g12.5)") arr(:,j,k)
       enddo
    enddo
  end subroutine
#:endfor
  
end module


