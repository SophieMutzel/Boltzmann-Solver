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
  integer(kind=ik)                     :: it, ierr, rank, nprocs, N,itask, L
  integer(kind=ik)                     :: io_error, istate, i, ier, neval, nit, nd
  real(kind=rk)                        :: rtol, atol, helper, T, abserr,s, H, epsabs, epsrel, rar, Tprime
  real(kind=rk), allocatable           :: Y(:), mas(:),aaxx(:,:)
  TYPE(VODE_OPTS)                      :: OPTIONS

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
  allocate(q(nrhs,params%N))
  L = int((params%z_max-params%z_start)/params%dz_plot)
  allocate(q_tot(nrhs+1,params%N,L+3))
  T = params%mx/10**params%z_start
  epsabs=1e-2_rk
  epsrel=1e-2_rk
  allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
  argsint%drhoa = params%drhoa
  call qags(rhop_over_rho,argsint,T,&
            T_RH, epsabs, epsrel, rar, abserr, neval, ier)
  if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
  deallocate(argsint%drhoa)
  Tprime = Ta(T,params,rar)
  s = ent(T,params)
  do i=1,params%N
    ! Y_chi=Yeq(T')
    q(1,i) = 0.0_rk
    q(2,i) = 0.0_rk!neq(T, params%ma, ga)
    q(3,i) = Tprime
  end do
  q_tot(1,:,1) = params%z_start
  q_tot(2,:,1) = q(1,:)/s
  q_tot(3,:,1) = q(2,:)/s

  dz = params%dz
  z  = params%z_start
  it = 0

  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "starting main time loop"
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
  conv_eps = 1e-3_rk
  z  = params%z_start
  do while ( z <= params%z_max .and. eps > conv_eps )!.or. z<=0.0_rk)
    it = it + 1
    zpdz = z + params%dz_plot
    CALL VODE_F90( region_freeze_in, nrhs*params%N, Y, z, zpdz, itask, &
                   istate, OPTIONS, params, argsint )
    q = reshape(Y,(/nrhs, params%N/))
    z = zpdz
    T = params%mx/10**z
    s = ent(T,params)
    H = Hub(T,rhoeq(T,params%ma,ga))
    q_tot(1,:,it) = z
    q_tot(2:3,:,it) = q(1:2,:)/s
    q_tot(4,:,it) = q(3,:)
    eps = abs((q_tot(1,1,it)-q_tot(1,1,it-1))/q_tot(1,1,it-1))

    write(*,*) z, eps
  end do

  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "main time loop over, writing results to files"
  end if

  ! write results to file
!  open (unit=97, file="temp/in_n.txt", status='old', action='write', position='append', iostat=io_error)
  do i=1,params%N
  !  call write_matrix("temp/"//exp2str(params%gaff(i))//exp2str(params%gaxx(i))//".txt",q_tot(:,i,:))
    call write_matrix("temp/"//trim(adjustl(params%file))//".txt",q_tot(:,i,1:it))
    call write_gnuplot(98, "temp/"//trim(adjustl(params%file)),Yxmxobs/params%mx, "fi")
!    if (io_error==0) then
!      write(97,*) params%gaxx(i)*params%gaff(i), params%gaxx(i), q_tot(2,i,it)
!    else
!      write(*,*) 'error', io_error,' while opening the file temp/in_n.txt'
!    end if
  end do
  !close(97)

  call MPI_Finalize(ierr)
end program main
