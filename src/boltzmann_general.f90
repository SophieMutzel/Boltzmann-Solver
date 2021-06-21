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
  integer(kind=ik)                     :: io_error, istate, i, ier, neval, nit, nd,nsqrt,j,k
  real(kind=rk)                        :: rtol, atol, T, abserr, Tprime, s, eps_old
  real(kind=rk), allocatable           :: Y(:), mas(:),aaxx(:,:), gaxx(:),kappa(:),gaxx_tot(:),gaff_tot(:)
  TYPE(VODE_OPTS)                      :: OPTIONS
  logical                              :: check

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
  call allocate_couplings( params )

  call ini_cons_to_params( params, argsint )

  argsint%mx = params%mx
  argsint%ma = params%ma

  argsint%g = 1.0_rk
  call choose_regime(params, argsint)
  if (params%regime .ne. "reannihilation") stop
  call initial_conditions( params, q, q_tot, rhs, argsint )

  dz = params%dz
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

  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL,H0=dz,HMAX=0.0001_rk,MXSTEP=100000)
  itask = 1
  istate = 1
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
  do while ( z <= params%z_max .and. eps > conv_eps .and. Tprime>1e-3_rk .or. z<=0.2_rk )
    !call rhs_contributions_general( nrhs*params%N, z, Y, params, argsint, rhs(:,:,it))
    it = it + 1
    zpdz = z + params%dz_plot
    CALL VODE_F90( general_rhs, nrhs*params%N, Y, z, zpdz, itask, &
                 istate, OPTIONS, params, argsint )
    q = reshape(Y,(/nrhs, params%N/))
    z = zpdz
    T = params%mx/10**z
    s = ent(T,params)
    q_tot(1,:,it) = z
    q_tot(2:nrhs,:,it) = q(1:2,:)/s
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

    !if (q_tot(nrhs+3,1,it)<q_tot(nrhs+3,1,it-1)) argsint%helper = .true.
    !eps = max(abs((q_tot(2,1,it)-q_tot(2,1,it-1))/q_tot(2,1,it-1)),abs((q_tot(3,1,it)-q_tot(3,1,it-1))/q_tot(3,1,it-1)))
    eps = abs((q_tot(2,1,it)-q_tot(2,1,it-1))/q_tot(2,1,it))

    write(*,*) z, eps, Tprime!q_tot(2,1,it), q_tot(4,1,it)
!    if (eps>eps_old .and. check) then
!      argsint%helper = .false.
!      check = .false.
!    end if
  end do

  !call check_BBN(q)
  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "main time loop over, writing results to files"
  end if

  ! write results to file
  open (unit=97, file="temp/reann_ratio5gaxx-2.txt", status='old', action='write', position='append', iostat=io_error)
  do i=1,params%N
  !  call write_matrix("temp/"//exp2str(params%gaff(i))//exp2str(params%gaxx(i))//".txt",q_tot(:,i,:))
    call write_matrix("temp/"//trim(adjustl(params%file))//".txt",q_tot(:,i,1:it))
    call write_matrix("temp/rhs_"//trim(adjustl(params%file))//".txt",rhs(:,i,1:it-1))
    call write_gnuplot(98, "temp/"//trim(adjustl(params%file)),Yxmxobs/params%mx, params%regime)
    call write_gnuplot_rhs(99, "temp/rhs_"//trim(adjustl(params%file)))
    if (io_error==0) then
      !write(97,*) params%gaxx(i)*params%gaff(i), params%gaxx(i), q_tot(2,i,it)
      write(97,*) params%gaff(i), params%ma, q_tot(2,i,it)
    else
      write(*,*) 'error', io_error,' while opening the file temp/in_n.txt'
    end if
  end do
  close(97)

  call MPI_Finalize(ierr)
end program main
