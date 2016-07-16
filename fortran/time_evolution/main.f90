! Ter Mai 10 10:15:50 BRT 2016
! author: gilberto medeiros nakamura
program time_evolution
  ! objective:
  !           given an initial condition, evaluates the
  !           probability vector (state vector) at time t
  !           for SIS model. A configuration is assembled
  !           using the N two-state configuration agents.
  ! procedure:
  !           a) contruct adjacency matrix
  !           b) construct time generator
  !             b.1) construct off-diagonal contributions
  !             b.2) construct diagonal contributions
  !           c) construct transition matrix
  !           d) apply the transition matrix on vector state
!  include 'mkl.f90'

  use markov_module
  implicit real*8(a-h,o-z)

!  include 'mkl.fi'

  integer            ::N,iconf,imax_size,isparse_size,i,j,k,lda,istatus
  integer,allocatable::A(:,:)
  integer,allocatable::ioccupation_vector(:),ioccupation_vector2(:)
  real*8 ,allocatable::H(:,:),prob(:),prob_new(:),mean_n(:),std_dev(:),rnorm(:)
  real*8 ,allocatable::values(:),prob_storage(:,:)
  integer,allocatable::icolumn_indx(:),irow_indx(:)
  real*8 ::params(0:1)
  integer::job(8)

  N=12
  imax_size=2**N


  params=0d0
  params(0)=1d0/(N*N)
  params(1)=(0.0d0)/N


  allocate(prob(0:imax_size-1))
  !initial condition
  prob=0d0  
  !prob=1d0/dble(imax_size) ! maximum entropy
  prob(1)=1d0 ! only agent 0 is infected

  
  allocate(A(0:N-1,0:N-1))
  allocate(H(0:imax_size-1,0:imax_size-1))

  print *,'generating transition matrix...'
  call create_meanfield_adjacency(N,A)  
  call create_time_generator(N,imax_size,params,A,H)
  ! H is now the transition matrix
  call create_transition_matrix(imax_size,H)
  print *,'   ... done'
  
  ! H is sparse ~ 10^-3 %
  ! first, compute the number of elements whose norm is
  ! greater than a minimum value 1d-16
  isparse_size=count(abs(H).gt.1d-16)
  print *,'sparticity :: ', real(isparse_size)/real(imax_size**2)

  ! values are the non-null values of H in CSR representation
  allocate(values(0:isparse_size-1),STAT=istatus)
  allocate(icolumn_indx(0:isparse_size-1),STAT=istatus)
  allocate(irow_indx(0:imax_size),STAT=istatus)


  ! use sparse blas level 2 to convert dense to sparse matrix
  job(1)=0 ! convert square to CSR
  job(2)=0 ! 0-index
  job(3)=0 ! 0-index
  job(4)=2 ! complete matrix
  job(5)=isparse_size
  job(6)=1
  lda=imax_size
  print *,'converting from dense to sparse representation...'
  call mkl_ddnscsr(job, imax_size, imax_size, H ,lda,values, icolumn_indx, irow_indx,istatus)
  !free memory
  if (allocated(H)) deallocate(H)
  print *,'   ... done    |    status =',istatus

  

  !perform time evolution
  print *,'starting time evolution...'
  imax_time = N**4
  allocate(prob_new(0:imax_size-1))
  allocate(rnorm(0:imax_time))
  allocate(mean_n(0:imax_time))
  allocate(std_dev(0:imax_time))
  allocate(ioccupation_vector(0:imax_size-1))
  allocate(ioccupation_vector2(0:imax_size-1))

  allocate(prob_storage(0:2,0:imax_time))

  prob_new=0d0
  mean_n =0d0
  std_dev=0d0
  call create_occupation_vector(imax_size,ioccupation_vector)
  ioccupation_vector2=ioccupation_vector**2d0
  mean_n(0) =dot_product(ioccupation_vector,prob)
  std_dev(0)=dot_product(ioccupation_vector2,prob)
  rnorm(0) = dot_product(prob,prob)

  prob_storage(0,0)=prob(0)
  prob_storage(1,0)=prob(5)
  prob_storage(2,0)=prob(imax_size-1)
  
  do iloop_time = 1,imax_time
     call mkl_cspblas_dcsrgemv('T',imax_size,values,irow_indx,icolumn_indx,prob,prob_new)
     prob = prob_new
     mean_n(iloop_time) =dot_product(ioccupation_vector,prob)
     std_dev(iloop_time)=dot_product(ioccupation_vector2,prob)
     rnorm(iloop_time)  =dot_product(prob,prob)


     prob_storage(0,iloop_time)=prob(0)
     prob_storage(1,iloop_time)=prob(5)
     prob_storage(2,iloop_time)=prob(imax_size-1)
  end do
  print *,'   ... done'

  std_dev=std_dev - (mean_n**2d0)
  std_dev=sqrt(abs(std_dev))/N
  mean_n=mean_n/N  
  do i=0,imax_time
     write(20,*) i,mean_n(i),std_dev(i)
     write(21,*) i,rnorm(i)
     write(22,*) prob_storage(:,i)
  end do

  

  
!!$  if (allocated(values)) deallocate(values,STAT=istatus)
!!$  print *,'   ... done',istatus
!!$ 
!!$  if (allocated(icolumn_indx)) deallocate(icolumn_indx,STAT=istatus)
!!$  print *,'   ... done',istatus
!!$  if (allocated(irow_indx)) deallocate(irow_indx,STAT=istatus)


  
end program time_evolution
