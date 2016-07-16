! Ter Mai 10 10:15:50 BRT 2016
! author: gilberto medeiros nakamura
module markov_module

contains
  subroutine create_meanfield_adjacency(N,A)
    !objective:
    !          create the adjacency matrix for fully connected graph
    !procedure:
    !          simply assign 1 to all entries but the main diagonal
    !variables:
    !          N = size
    !          A = adjacency matrix ( i=0..N-1, j=0..N-1)
    integer::N
    integer::A(0:N-1,0:N-1)
    integer::i
    A = 1
    do i=0,N-1
       A(i,i)=0
    end do
  end subroutine create_meanfield_adjacency
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine create_random_adjacency(N,p,A)
    !objective:
    !          create the adjacency matrix for random graph (uniform
    !          distribution)
    !procedure:
    !          bonds appear with probability p
    !variables:
    !          N = size
    !          A = adjacency matrix ( i=0..N-1, j=0..N-1)
    !          p = connection probability
    integer::N
    integer::A(0:N-1,0:N-1)
    real*8 ::p,B(0:N-1,0:N-1)
    integer::i,j
    call random_number(B)
    where (B.lt.p)
       A=1
    end where
    do j=0,N-1
       do i=0,N-1
          A(i,j)=A(j,i)
       end do
    end do
    do i=0,N-1
       A(i,i)=0
    end do
    
  end subroutine create_random_adjacency
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine create_time_generator(N,imax_size,params,A,H)
    implicit real*8(a-h,o-z)
    !objective:
    !          create the time generator H, given the adjacency matrix A (NxX)
    !          and the parameters params for SIS model.
    !procedure:
    !          H is written as
    !          H = params(0)*sum_{j k}A_{j k}{ n_k *(1-n_j) - \sigma_k^+ n_j}+
    !            + params(1)*sum_{k}{n_k - \sigma_k^-}
    !          Operators act on configuration states, whose representation are
    !          binary numbers \mu = n_0*2^0 +n_1*2^1 +n_2*2^2+...+n_{N-1}*2^{N-1}.
    !          The action \sigma_{k}^{+-} over states simply change the bit.
    !variables:
    !          N        = size
    !          imax_size= 2**N -1
    !          A        = adjacency matrix ( i=0..N-1, j=0..N-1)
    !          params   = parameters. Here, params(0) = transmission rate and
    !                     params(1) = cure rate
    !          H        = time generators
    integer::N,imax_size
    integer::A(0:N-1,0:N-1)
    real*8::params(0:1) 
    real*8::H(0:imax_size-1,0:imax_size-1)
    real*8::rtmp
    integer::iconf,inew_conf,k,j,i,ioccupation,idelta,ioccupation_j,ioccupation_i
    !comment:
    !        one viable implementation sorts configurations by occupation number.
    !        The approach improves access time and allows better OMP parallelization.

    !>>first: off-diagonal contributions<<
    
    !loop over configuration iconf.
    do iconf=0,imax_size-1
       !apply the off-diagonal \sigma^-_k over each site k in iconf
       !store new configuration in inew_conf.
       !obs: access time for H could be improved if inew_conf is stored in
       !     a list, which would be sorted and then read.
       !     However, since each local off-diagonal operator produce at most
       !     one new configurations, the overhead and posterior sorting would
       !     actually increase time effort.
       do k=0,N-1
          !apply_cure produce the action of cure operators over iconf
          !resulting in inew_conf. If the action produces null vector,
          !then idelta = 0
          call apply_cure(k,iconf,inew_conf,idelta)
          H(inew_conf,iconf) = H(inew_conf,iconf) - idelta * params(1) 
       end do
    end do

    !loop over configuration iconf.
    do iconf=0,imax_size-1
       !apply the off-diagonal \sigma^+_k n_j over each site k in iconf
       !store new configuration in inew_conf.
       !obs: access time for H could be improved if inew_conf is stored in
       !     a list, which would be sorted and then read.
       !     However, since each local off-diagonal operator produce at most
       !     one new configurations, the overhead and posterior sorting would
       !     actually increase time effort.
       do j=0,N-1
          ! check if the j-th agent is infected
          !ioccupation = ibits(iconf,j,1)
          ioccupation = btest(iconf,j)
          if (ioccupation) then
             do k=0,N-1
                ! apply_infection works similarly to apply_cure
                ! the only difference being the \sigma^+_k operator
                call apply_infection(k,iconf,inew_conf,idelta)             
                H(inew_conf,iconf) = H(inew_conf,iconf) - idelta * params(0) * A(k,j)
             end do
          end if
       end do
    end do

    !diagonal contribution
    do iconf=0,imax_size-1
       rtmp=0d0
       do j=0,N-1
          ioccupation_j=btest(iconf,j)
          if (ioccupation_j) then
             !add the cure contribution
             rtmp = rtmp + params(1)
             do i=0,N-1
                ioccupation_i = btest(iconf,i)
                if (.not.ioccupation_i) then
                   rtmp = rtmp + params(0)*A(i,j)
                end if
             end do
          end if
       end do
       H(iconf,iconf)=rtmp
    end do   
    
  end subroutine create_time_generator
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine apply_cure(k,iconf,inew_conf,idelta)
    !objective:
    !          apply the \sigma_k^- operator on iconf producing inew_conf
    !procedure:
    !          set k-th bit to null = inew_conf.
    !          if inew_conf XOR iconf = 0 then return < 0 (null vector)
    !variables:
    !          k = position/label
    !          iconf    = configuration  (input) (column)
    !          inew_conf= new configuration (output)
    !          idelta   = change flag indicator
    integer::k,iconf,inew_conf,idelta
    inew_conf = ibclr(iconf,k)
    !evaluate bit-difference
    !if idelta = 0, inew_conf = iconf or iconf = 0
    !if idelta = 1, inew_conf <> iconf
    idelta = min(ieor(inew_conf,iconf), 1)
  end subroutine apply_cure
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine apply_infection(k,iconf,inew_conf,idelta)
    !objective:
    !          apply the \sigma_k^+ operator on iconf producing inew_conf
    !procedure:
    !          set k-th bit to null = inew_conf.
    !          if inew_conf XOR iconf = 0 then return < 0 (null vector)
    !variables:
    !          k = position/label
    !          iconf    = configuration  (input) (column)
    !          inew_conf= new configuration (output)
    !          idelta   = change flag indicator
    integer::k,iconf,inew_conf,idelta
    inew_conf = ibset(iconf,k)
    !evaluate bit-difference
    !if idelta = 0, inew_conf = iconf or iconf = 0
    !if idelta = 1, inew_conf <> iconf
    idelta = min(ieor(inew_conf,iconf), 1)
  end subroutine apply_infection
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine create_transition_matrix(imax_size,H)
    integer::imax_size
    real*8 ::H(0:imax_size-1,0:imax_size-1)
    integer::iconf
    H=-H  ! multiply by time interval
    do iconf=0,imax_size-1
       H(iconf,iconf)=1d0 + H(iconf,iconf)
    end do
  end subroutine create_transition_matrix
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine create_occupation_vector(imax_size,n_vector)
    !objective:
    !          the idea is to use the diagonal elements to represent
    !          the occupation operator (global)
    !procedure:
    !          use popcnt(iconf) over configuration iconf to obtain
    !          the number of infected agents
    !variables:
    !          imax_size = hilbert space dimension
    !          n_vector  = output (diagonal of n operator)
    integer:: imax_size
    integer:: n_vector(0:imax_size-1)
    integer:: iconf
    do iconf=0, imax_size-1
       n_vector(iconf)=popcnt(iconf)
    end do    
  end subroutine create_occupation_vector
end module markov_module
