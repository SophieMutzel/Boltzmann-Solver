program main
!------------------------------------------------------------------------------!
  use mpi
  use module_params
  use module_utils
  use module_cosmo
  use module_xsecs
  use module_rhs
  use DVODE_F90_M
!------------------------------------------------------------------------------!
  implicit none
  ! global parameters
  real(kind=rk)                        :: z, dz, eps, zpdz, conv_eps
  type (type_params)                   :: params
  type (type_argsint)                  :: argsint
  real(kind=rk), allocatable           :: q(:,:)
  real(kind=rk), allocatable           :: q_tot(:,:,:), rhs(:,:,:)
  integer(kind=ik)                     :: it, ierr, rank, nprocs, N,itask
  integer(kind=ik)                     :: io_error, istate, i, ier, neval, nit, nd,j,k
  real(kind=rk)                        :: rtol, atol, T, abserr, Tprime, s, eps_old,n0, test
  real(kind=rk)                        :: gam_agff, ffa, gam_xxff, sv_aaxx, e, ma, mx
  real(kind=rk), allocatable           :: Y(:), gaxx(:), gaff(:)
  TYPE(VODE_OPTS)                      :: OPTIONS
  logical                              :: check
  real(kind=rk), dimension(100) :: stest
!------------------------------------------------------------------------------!
  ! init mpi
  call MPI_Init(ierr)
  ! process rank
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  params%rank = rank
  ! process number
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  params%nprocs = nprocs
  ! output MPI status
  params%BOLTZMANN_COMM = MPI_COMM_WORLD
  call set_mpi_comm_global(MPI_COMM_WORLD)

  if (rank==0) then
    write(*,'(80("_"))')
    write(*, '("number of processes:", i5)') params%nprocs
  end if
  nrhs = 3
  call get_params( params )
  ! create parameter struct
  call ini_cons_to_params( params, argsint )
  argsint%mx = params%mx
  argsint%ma = params%ma
  ! if we want to build a grid of couplings...
  if (params%grid) then
    params%N = 1_ik
    allocate(params%gaff(params%N))
    allocate(params%gaxx(params%N))
    call build_grid(params, argsint, gaff, gaxx)
    open (unit=97, file="temp/gaffgaxx.txt", status='old', action='write', position='append', iostat=io_error)
    do j=1,size(gaff)
      if (io_error==0) then
        write(97,100) params%ma, params%mx, log10(gaff(j)), log10(gaxx(j))
        100 format( F10.3 , F10.3, F10.2, F10.2)
      else
        write(*,*) 'error', io_error,' while opening the file temp/gaffgaxx.txt'
      end if
    end do
    close(97)
    ! call freeze_in_grid(params, argsint)
  else
    call allocate_couplings( params )
  end if
  argsint%g = 1.0_rk
  ! which production regime are we in (freeze-in, reannihilation, sequential freeze-in, freeze-out)?
  call choose_regime(params, argsint)
  ! calculate intial conditions
  call initial_conditions( params, q, q_tot, rhs, argsint )
  ! calculate af->gf in advance and later interpolate to speed up computation
  call aftogf( params, argsint)

  ! time step
  dz = params%dz
  ! initial time
  z  = params%z_start
  it = 0
  if (rank==0) then
    write(*,'(80("-"))')
    write(*,*) "starting main time loop"
    write(*,*) "Regime: ", params%regime
    write(*,'(80("-"))')
  end if

  atol = 1e-20_rk ! absolute tolerance
  rtol = 1e-10_rk ! relative tolerance

  ! Set options for dvode
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL,H0=dz,HMAX=0.001_rk,MXSTEP=100000)
  itask = 1
  istate = 1
  ! array for final relic densities/HS temperature
  allocate(Y(nrhs*params%N))
  Y = reshape(q, (/nrhs*params%N/) )
  it = 1

  eps = 1.0_rk
  conv_eps = 5e-4_rk
  z  = params%z_start
  argsint%helper = .false.
  Tprime = 1.0_rk
  check = .true.
  eps_old=1.0_rk
  ! main time loop
  do while ( z<=params%z_max .and. eps > conv_eps  .or. Tprime>params%mx/40.0_rk )
    ! routine to plot the reaction rates of the rhs
    call rhs_contributions_general( nrhs*params%N, z, Y, params, argsint, rhs(:,:,it))
    it = it + 1
    zpdz = z + params%dz_plot
    ! if we are in reannihilation or freeze-out, solve for log_10(n) for better stability
    if (params%regime == "reannihilation" .or. params%regime == "freeze-out") then
      CALL VODE_F90( boltzmann_logn, nrhs*params%N, Y, z, zpdz, itask, &
                 istate, OPTIONS, params, argsint )
    else
      ! for freeze-in solve for n since initial condition is 0
      Y(1:2) = 10**Y(1:2)
      if (it==1) Y(1:2) = 0.0_rk
      CALL VODE_F90( general_rhs, nrhs*params%N, Y, z, zpdz, itask, &
                 istate, OPTIONS, params, argsint )
      Y(1:2) = log10(Y(1:2))
    end if
    q = reshape(Y,(/nrhs, params%N/))
    ! new value of z
    z = zpdz
    ! temperature of current step
    T = params%mx/10**z
    ! entropy of current step (use Y:=n/T^3 to get rid of dof)
    s = T*T*T!ent(T,params)
    ! save result of current time step in q_tot
    q_tot(1,:,it) = z
    !q_tot(2:nrhs,:,it) = q(1:2,:)/s
    q_tot(2:nrhs,:,it) = q(1:2,:)-log10(s)
    select case (params%regime)
    case ("reannihilation")
      Tprime = q(3,1)
    case ("freeze-in")
      Tprime = T
    case ("seq-freeze-in")
      Tprime = q(3,1)
    case ("freeze-out")
      Tprime = T
    end select
    q_tot(nrhs+1,:,it) = neq(T, params%mx, gDM)/s
    do i=1,params%N
      q_tot(nrhs+2,i,it) = neq(Tprime, params%mx, gDM)/s
      q_tot(nrhs+3,i,it) = neq(Tprime, params%ma, ga)/s
    end do
    q_tot(nrhs+4,:,it) = Tprime

    ! convergence criterion
    eps = abs((q_tot(2,1,it)-q_tot(2,1,it-1))/q_tot(2,1,it))
    write(*,*) z, eps, Tprime!q_tot(2,1,it), q_tot(4,1,it)
  end do
  !call check_BBN(q,T,params)
  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "main time loop over, writing results to files"
  end if
  q_tot(nrhs+1:nrhs+3,:,1:it)=log10(q_tot(nrhs+1:nrhs+3,:,1:it))
  ! DM density today
  n0 = omegax*rhocrit/params%mx
  !write(*,*) 10.0**q(1,1)/ent(T,params)*s0*params%mx/rhocrit, log10(n0/s0*ent(T,params)/s/2.0_rk)
  !write(*,*) q_tot(2,1,it), log10(n0/s0*ent(T,params)/s/2.0_rk)
  ! write results to file
  open (unit=97, file="temp/output.txt", status='old', action='write', position='append', iostat=io_error)
  do i=1,params%N
    call write_matrix("temp/"//trim(adjustl(params%file))//".txt",q_tot(:,i,1:it))
    call write_matrix("temp/rhs_"//trim(adjustl(params%file))//".txt",rhs(:,i,1:it-1))
    call write_gnuplot(98, "temp/"//trim(adjustl(params%file)),log10(n0/s0*ent(T,params)/s/2.0_rk), params%regime)
    call write_gnuplot_rhs(99, "temp/rhs_"//trim(adjustl(params%file)))
    if (io_error==0) then
      !write(97,*) params%gaxx(i)*params%gaff(i), params%gaxx(i), q_tot(2,i,it)
      !write(97,*) params%ma, params%gaff(i), params%gaxx(i), 10.0**q(1,i)/ent(T,params)*s0*params%mx/rhocrit!q_tot(2,i,it)
      write(97,*) params%gaff(i)*params%gaxx(i), params%gaxx(i), 10.0**q(1,i)/ent(T,params)*s0*params%mx/rhocrit
      !write(97,*) 10.0**q(2,i)/ent(T,params), params%ma, params%gaff(i), params%gaxx(i)
    else
      write(*,*) 'error', io_error,' while opening the file temp/output.txt'
    end if
  end do
  close(97)
  !14 continue
  call MPI_Finalize(ierr)
end program main
