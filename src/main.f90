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
  integer(kind=ik)                     :: io_error, istate, i, ier, neval, nit, nd
  real(kind=rk)                        :: rtol, atol, helper, T, test, abserr, Tprime, s, ztest, Told
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

  argsint%g = 1.0_rk

  allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
  argsint%drhoa = params%drhoa
  params%regime="reannihilation"
  dz = params%dz
  z  = params%z_start
  it = 0
  ztest = 0.0_rk
  Told = T_RH
!  do while ( z <= params%z_max)
!!    !call test_bezier(params%geff_HS(1,:),params%geff_HS(2,:),z,test,nd-1,params%A,params%B)
!    zpdz = z +params%dz_plot
!    T = params%mx/10**z
!    call qags(rhop_over_rho,argsint,T_RH,&
!              T, 1e-30_rk, 1e-30_rk, test, abserr, neval, ier)
!    !Told = T
!    !ztest = test + ztest
!!    test = geff_rho(T)
!!    call geffSM(T,params,Tprime)
!    !Tprime = Ta(T,params)
!    ! HS interaction
!    !argsint%g = params%gaxx(1)
!    !call sigmav( Tprime, params, argsint, "aaxx", test )
!!    !test = 2.46743 - 0.900703 *z - 0.426853 *z**2 + 0.344933 *z**3 + 0.241269 *z**4 - 1.64352 *tanh(2.49447*z)
!!    !call interp_linear(size(x_eval), x_eval,y_eval,z,test)!log10(Tprime), test)
!    z = zpdz
!    !s = ent(T,params)
!
!    !write(*,*) z, test, Hub( T, params ), ent( T, params )
!    write(*,*) z, T, test
!  end do
!    deallocate(argsint%drhoa)
!  z = params%z_start
!stop

  call initial_conditions( params, q, q_tot, rhs, argsint )

  ! check if equilibrium in HS attained
  call sigmav( q(3,1), params, argsint, "aaxx", test )
  if (test*params%gaxx(1)*params%gaxx(1)*params%gaxx(1)&
      *params%gaxx(1)*neq(q(3,1), params%ma, ga)<Hub(params%mx/10**params%z_start,0.0_rk ))&
      write(*,*) "Eq not attained!"

  dz = params%dz
  z  = params%z_start
  it = 0
  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "starting main time loop"
  end if

  atol = 1e-20_rk ! absolute tolerance
  rtol = 1e-10_rk ! relative tolerance

  !options = set_opts( DENSE_J = .true.,USER_SUPPLIED_JACOBIAN= .false.,&
                    !RELERR=rtol,ABSERR=atol )
                    !METHOD_FLAG=25,LOWER_BANDWIDTH=nrhs-1, UPPER_BANDWIDTH=nrhs-1,RELERR=rtol,ABSERR=atol )
  !options = set_opts( METHOD_FLAG=25, BANDED_J = .true.,USER_SUPPLIED_JACOBIAN= .false.,&
  !                    RELERR=rtol,ABSERR=atol, LOWER_BANDWIDTH=nrhs-1, UPPER_BANDWIDTH=nrhs-1 )
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL,H0=dz,HMAX=0.0001_rk,MXSTEP=100000)
  itask = 1
  istate = 1
  allocate(Y(nrhs*params%N))
  Y = reshape(q, (/nrhs*params%N/) )
  it = 1

  eps = 1.0_rk
  conv_eps = 1e-3_rk
  z  = params%z_start
  argsint%helper = .false.
  do while ( z <= params%z_max .and. eps > conv_eps .or. z<=0.0_rk)
    call rhs_contributions_in_n( nrhs*params%N, z, Y, params, argsint, rhs(:,:,it), "reann" )
    it = it + 1
    zpdz = z + params%dz_plot

    CALL VODE_F90( region3a_in_n, nrhs*params%N, Y, z, zpdz, itask, &
                   istate, OPTIONS, params, argsint )
    q = reshape(Y,(/nrhs, params%N/))
    z = zpdz
    T = params%mx/10**z
    s = ent(T,params)
    q_tot(1,:,it) = z
    q_tot(2:nrhs,:,it) = q(1:2,:)/s
    Tprime = q(3,1)
    q_tot(nrhs+1,:,it) = neq(T, params%mx, gDM)/s
    do i=1,params%N
      q_tot(nrhs+2,i,it) = neq(Tprime, params%mx, gDM)/s
      q_tot(nrhs+3,i,it) = neq(Tprime, params%ma, ga)/s
    end do
    q_tot(nrhs+4,:,it) = Tprime!params%mx/(sqrt(params%gaff)*Tprime)

    if (q_tot(nrhs+3,1,it)<q_tot(nrhs+3,1,it-1)) argsint%helper = .true.
    eps = max(abs((q_tot(2,1,it)-q_tot(2,1,it-1))/q_tot(2,1,it-1)),abs((q_tot(3,1,it)-q_tot(3,1,it-1))/q_tot(3,1,it-1)))

    !write(*,*) z, eps, Tprime, T
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
    call write_matrix("temp/rhs_"//trim(adjustl(params%file))//".txt",rhs(:,i,1:it-1))
    call write_gnuplot(98, "temp/"//trim(adjustl(params%file)),Yxmxobs/params%mx, "reannihilation")
    call write_gnuplot_rhs(99, "temp/rhs_"//trim(adjustl(params%file)))
!    if (io_error==0) then
!      write(97,*) params%gaxx(i)*params%gaff(i), params%gaxx(i), q_tot(2,i,it)
!    else
!      write(*,*) 'error', io_error,' while opening the file temp/in_n.txt'
!    end if
  end do
!  close(97)

  call MPI_Finalize(ierr)
end program main
