
subroutine FCLS_single_pixel(rr1, AA, abundance)
  use dm
  use dm_op
  use dm_type

  implicit none
  real(kind=8), pointer :: rr1(:)
  real(kind=8), pointer :: AA(:, :)
  real(kind=8), pointer :: abundance(:)
  real(kind=8), allocatable :: AA1D(:), lagra_array(:), &
       tmp_array1(:), tmp_array2(:)
  integer :: ierr
  integer :: i, u, v, m, n, numloop
  integer :: size_r1(1), size_A(2)

  type(Matrix) :: r1, A, sample, EtE, eA, E, er, f, Etf
  type(Matrix) :: ls, x2, x2old, Z, L, L1
  type(Matrix) :: weights, xerow
  type(Matrix) :: iEtEOne, iEtE, One, lagra, maxneg  
  
  real(kind=8) :: max_val
  integer :: max_pos(3), iter, iter_on_whole
  real(kind=8) :: ee, lamdiv2(1)
  real(kind=8) :: sumiEtEOne, tol
  integer, allocatable :: zz(:)
  integer :: lastrow, lastcol, ab, iz
  real(kind=8) :: lagra_lastrow(1)
  size_r1 = shape(rr1)
  size_A  = shape(AA)



  r1 = dm_zeros(size_r1(1), 1, 1)
  A  = dm_zeros(size_A(1), size_A(2), 1)
  maxneg = dm_zeros(1, 1, 1)
  
  allocate(AA1D(A%nx*A%ny))

  call dm_setvalues(r1, (/(i,i=0,r1%nx-1)/),(/0/),(/0/), rr1, ierr)
  AA1D = reshape(AA, shape(AA1D))
  
  call dm_setvalues(A, (/(i,i=0,A%nx-1)/),(/(i,i=0,A%ny-1)/), (/0/), AA1D, ierr)

  !call dm_view(A, ierr)
  !call dm_view(r1, ierr)

  ! numloop=size(A,2);
  numloop = A%ny

  ! ee=1/10/max(r1);
  call dm_max(r1, max_val, max_pos, ierr)
  ee = 1.0 / 10. / max_val

  ! eA=e*A;
  eA = ee * A

  ! E=[ones(1,numloop);eA];
  E = dm_ones(1, numloop, 1) .xj. eA

  ! EtE=E'*E;
  EtE = dm_trans(E) * E
  ! [m,n] = size(EtE);
  m = EtE%nx
  n = EtE%ny

  !call dm_view(EtE, ierr)

  ! One=ones(m,1);
  One = dm_ones(m, 1, 1)

  ! iEtE=inv(EtE);
  iEtE = dm_inverse(EtE)

  !call dm_view(iEtE * EtE, ierr)

  ! iEtEOne=iEtE*One;
  iEtEOne = iEtE * One

  ! sumiEtEOne=sum(iEtEOne);
  sumiEtEOne = dm_sum_all(iEtEOne)

  !   ! weights=diag(iEtE);
  weights = dm_getdiag(iEtE)

  ! c=0;
  !c = 0

  ! sample=r1;
  sample = r1

  ! er=e*sample;
  er = ee * sample

  ! f=[1;er];
  f = dm_ones(1, 1, 1) .xj. er

  ! Etf=E'*f;
  Etf = dm_trans(E) * f

  tol=1e-7;

  call dm_view(Etf, ierr)

  ! %fcls1a

  ! %%%% THIS IS lamdiv2
  ! ls=iEtE*Etf;
  ls = iEtE * Etf

  ! ! lamdiv2=-(1-(ls'*One))/sumiEtEOne;
  !call dm_view(-(1-(dm_trans(ls)*One))  * (1.0/sumiEtEOne), ierr)
  call dm_getvalues(-(1-(dm_trans(ls)*One))  * (1.0/sumiEtEOne), &
       (/0/), (/0/), (/0/), lamdiv2, ierr)

  ! x2=ls-lamdiv2*iEtEOne;
  x2 = ls - lamdiv2(1) * iEtEone

  ! x2old=x2;
  x2old = x2

  ! iter_on_whole=0;
  iter_on_whole = 0

  call dm_view(x2, ierr)

  ! if (any(x2<-tol))&&(iter_on_whole<=10*m)
  if(dm_sum_all(x2 < -tol) > 0 .and. iter_on_whole <= 10*m) then

     !iter_on_whole=iter_on_whole+1;
     iter_on_whole = iter_on_whole + 1

     !Z=zeros(m,1);
     Z = dm_zeros(m, 1, 1)

     !iter=0;
     iter = 0

     !  while(any(x2<-tol) && iter <(m))
     do while((dm_sum_all(x2 < -tol)>0.1) .and. (iter < m))
        ! Z(x2<-tol)=1;
        Z = (Z .em. (1-(x2<-tol))) + (x2<-tol)
        !call dm_view(x2 < -tol, ierr)
        !call dm_view(Z, ierr)

        !zz=find(Z);
        zz = dm_find2(Z)

        ! x2=x2old;              % Reset x2
        x2 = x2old

        ! L=iEtE(zz,zz);
        L = dm_getsub(iEtE, zz, zz, (/0/))

        ! ab=size(zz);
        ab = size(zz)

        ! lastrow=ab(1)+1;
        lastrow = ab + 1

        !lastcol=lastrow;
        lastcol = lastrow

        ! L(lastrow,1:ab(1))=(iEtE(:,zz)'*One)';      
        allocate(tmp_array2(ab), tmp_array1(ab*ab))

        L1 = dm_zeros(ab+1, ab+1, 1)

        call dm_getvalues(L,  (/(i,i=0,ab-1)/), (/(i,i=0,ab-1)/), (/0/), &
             tmp_array1, ierr)
        call dm_setvalues(L1, (/(i,i=0,ab-1)/), (/(i,i=0,ab-1)/), (/0/), &
             tmp_array1, ierr)
        call dm_getvalues(dm_trans(dm_getsub(iEtE, (/(i,i=0,iEtE%nx-1)/), &
             zz, (/0/)) ) * One, (/(i,i=0,ab-1)/), (/0/), (/0/), &
             tmp_array2, ierr)
        call dm_setvalues(L1, (/lastrow-1/), (/(i,i=0,ab-1)/), (/0/), &
             tmp_array2, ierr)

        L = L1
        deallocate(tmp_array1, tmp_array2)

        !L(1:ab(1),lastcol)=iEtEOne(zz);
        !L(lastrow,lastcol)=sumiEtEOne;      
        allocate(tmp_array1(size(zz)))
        call dm_getvalues(iEtEOne, (/(i,i=0,ab-1)/), (/0/), (/0/), &
             tmp_array1, ierr)
        call dm_setvalues(L, (/(i,i=0,ab-1)/), (/lastcol-1/), (/0/), &
             tmp_array1, ierr)
        call dm_setvalues(L, (/lastrow-1/), (/lastcol-1/), (/0/), &
             (/sumiEtEOne/),ierr)
        deallocate(tmp_array1)

        !call dm_view(L, ierr)

        !  xerow=x2(zz);
        !  xerow(lastrow,1)=0;
        ! NOTE : the size of xerow changes!!
        xerow = dm_zeros(lastrow, 1, 1)
        allocate(tmp_array1(size(zz)))
        call dm_getvalues(x2, zz, (/0/), (/0/), tmp_array1, ierr)
        call dm_setvalues(xerow, (/(i,i=0,size(zz)-1)/), (/0/), (/0/), &
             tmp_array1, ierr)
        call dm_setvalues(xerow, (/lastrow-1/), (/0/), (/0/), (/0./), ierr)
        !call dm_view(xerow, ierr)
        deallocate(tmp_array1)

        ! lagra=L\xerow;
        lagra = dm_solve(L, xerow)

        if(allocated(lagra_array)) deallocate(lagra_array)
        allocate(lagra_array(ab))
        lagra_array = 0
        call dm_getvalues(lagra, (/(i,i=0,ab-1)/), (/0/), (/0/),&
             lagra_array, ierr)

        ! print*, "xerow="
        ! call dm_view(xerow, ierr)
        ! print*, "lagra="
        ! call dm_view(lagra, ierr)
        
        !  while (any(lagra(1:ab(1))>0))   % Reset lagrange multipliers
        do while(any(lagra_array > 0))
           ! maxneg=weights(zz).*lagra(1:ab(1));
           maxneg = dm_getsub(weights, zz, (/0/), (/0/)) .em. &
                dm_getsub(lagra, (/(i,i=0,ab-1)/), (/0/), (/0/))

           ! [yz,iz]=max(maxneg);   % Remove the most positive
           call dm_max(maxneg, max_val, max_pos, ierr)
           iz = max_pos(1)

           ! Z(zz(iz))=0;
           call dm_setvalues(Z, (/zz(iz)/), (/0/), (/0/), (/0/), ierr)

           ! zz=find(Z);           % Will always be at least one (prove)
           deallocate(zz)
           zz = dm_find2(Z)
           ! L=iEtE(zz,zz);
           L = dm_getsub(iEtE, zz, zz, (/0/))

           ! ab=size(zz);
           ab = size(zz)

           ! lastrow=ab(1)+1;
           lastrow = ab + 1

           ! lastcol=lastrow;
           lastcol = lastrow

           ! L(lastrow,1:ab(1))=(iEtE(:,zz)'*One)';
           allocate(tmp_array2(ab), tmp_array1(ab*ab))

           L1 = dm_zeros(ab+1, ab+1, 1)

           call dm_getvalues(L,  (/(i,i=0,ab-1)/), (/(i,i=0,ab-1)/), (/0/), &
                tmp_array1, ierr)
           call dm_setvalues(L1, (/(i,i=0,ab-1)/), (/(i,i=0,ab-1)/), (/0/), &
                tmp_array1, ierr)
           call dm_getvalues(dm_trans(dm_getsub(iEtE, (/(i,i=0,iEtE%nx-1)/), &
                zz, (/0/)) ) * One, (/(i,i=0,ab-1)/), (/0/), (/0/), &
                tmp_array2, ierr)
           call dm_setvalues(L1, (/lastrow-1/), (/(i,i=0,ab-1)/), (/0/), &
                tmp_array2, ierr)

           L = L1
           deallocate(tmp_array1, tmp_array2)

           ! L(1:ab(1),lastcol)=iEtEOne(zz);
           ! L(lastrow,lastcol)=sumiEtEOne;
           allocate(tmp_array1(size(zz)))
           call dm_getvalues(iEtEOne, (/(i,i=0,ab-1)/), (/0/), (/0/), &
                tmp_array1, ierr)
           call dm_setvalues(L, (/(i,i=0,ab-1)/), (/lastcol-1/), (/0/), &
                tmp_array1, ierr)
           call dm_setvalues(L, (/lastrow-1/), (/lastcol-1/), (/0/), &
                (/sumiEtEOne/),ierr)

           deallocate(tmp_array1)

           ! xerow=x2(zz);
           ! xerow(lastrow,1)=0;
           xerow = dm_zeros(lastrow, 1, 1)
           allocate(tmp_array1(size(zz)))
           call dm_getvalues(x2, zz, (/0/), (/0/), tmp_array1, ierr)

           call dm_setvalues(xerow, (/(i,i=0,size(zz)-1)/), (/0/), (/0/), &
                tmp_array1, ierr)

           call dm_setvalues(xerow, (/lastrow-1/), (/0/), (/0/), (/0./), ierr)
           !call dm_view(xerow, ierr)
           deallocate(tmp_array1)

           ! lagra=L\xerow;
           lagra = dm_solve(L, xerow)
           lagra_array = 0
           call dm_getvalues(lagra, (/(i,i=0,ab-1)/), (/0/), (/0/), &
                lagra_array, ierr)
        enddo

        ! %problem with lamscls zz may be null
        ! if ~isempty(zz)
        ! x2=x2-iEtE(:,zz)*lagra(1:ab(1))-lagra(lastrow)*iEtEOne;
        if(size(zz)/=0) then
           call dm_getvalues(lagra, (/lastrow-1/), (/0/), (/0/), &
                lagra_lastrow, ierr)
           x2 = x2 - (dm_getsub(iEtE, (/(i,i=0,iEtE%nx-1)/), zz, (/0/))* &
                dm_getsub(lagra, (/(i,i=0,ab-1)/), (/0/), (/0/))) -&
                lagra_lastrow(1) * iEtEOne;
        endif
        iter=iter+1;
        deallocate(zz)      
     enddo
  endif

  if(allocated(lagra_array)) deallocate(lagra_array)
  deallocate(AA1D)

  !abundance=x2;  
  call dm_getvalues(x2, (/(i,i=0,x2%nx-1)/), (/0/), (/0/), abundance, ierr)

  call dm_destroy(r1, ierr)
  call dm_destroy(A, ierr)
  call dm_destroy(sample, ierr)
  call dm_destroy(EtE, ierr)
  call dm_destroy(eA, ierr)
  call dm_destroy(E, ierr)
  call dm_destroy(er, ierr)
  call dm_destroy(f, ierr)
  call dm_destroy(Etf, ierr)
  call dm_destroy(ls, ierr)
  call dm_destroy(x2, ierr)
  call dm_destroy(x2old, ierr)
  call dm_destroy(Z, ierr)
  call dm_destroy(L, ierr)
  call dm_destroy(L1, ierr)
  call dm_destroy(weights, ierr)
  call dm_destroy(xerow, ierr)
  call dm_destroy(iEtEOne, ierr)
  call dm_destroy(iEtE, ierr)
  call dm_destroy(One, ierr)
  call dm_destroy(lagra, ierr)
  call dm_destroy(maxneg, ierr)

end subroutine FCLS_single_pixel

program main
  use dm
  use dm_op
  use dm_type
  implicit none
  include 'mpif.h'
  interface
     subroutine   FCLS_single_pixel(rr1, AA1, abundance1)
       real(kind=8), pointer :: rr1(:)
       real(kind=8), pointer :: AA1(:, :)
       real(kind=8), pointer :: abundance1(:)
     end subroutine FCLS_single_pixel
  end interface

  real(kind=8), pointer :: r1(:)
  real(kind=8), pointer :: A(:, :)
  real(kind=8), pointer :: abundance(:)
  real(kind=8), allocatable :: tt(:)
  integer :: N, i, ierr
  integer :: rank, size

  call dm_init1(ierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  
  N = 20
  allocate(r1(N), A(N, N), abundance(N))

  do i = 1,N
     r1(i) = i
     A(i, i) = 2.0;
  enddo

  call FCLS_single_pixel(r1, A, abundance)

  OPEN(UNIT=12, FILE="abundance.txt", ACTION="write", STATUS="replace")
  if(rank == 0) then
     !print*, "abundance="
     !write(*, "(T4, ES15.8)") abundance
     write(12, "(T1, ES15.8)") abundance
  endif

  deallocate(r1, A, abundance)
  call dm_finalize(ierr)  
end program main

