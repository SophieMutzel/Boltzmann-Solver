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
  real(kind=rk)                        :: z, dz, eps, zpdz
  type (type_params)                   :: params
  type (type_argsint)                  :: argsint
  real(kind=rk), allocatable           :: q(:,:), q_new(:,:)
  real(kind=rk), allocatable           :: q_tot(:,:,:)
  integer(kind=ik)                     :: it, ierr, rank, nprocs, N,itask, istate, i
  real(kind=rk)                        :: rtol, atol, helper, T
  real(kind=rk), allocatable           :: Y(:)
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

  ! create parameter struct
  call allocate_couplings( params)

  call ini_cons_to_params( params )

  call initial_conditions( params, q, q_new, q_tot )

  argsint%mx = params%mx
  argsint%ma = params%ma

  dz = params%dz
  z  = params%z_start
  it = 0

!  if (rank==0) then
!    write(*,'(80("_"))')
!    write(*,*) "optimizing dz"
!  end if
!  N = size(q)
!  eps = 1d-3
  ! -------- MAIN TIME LOOP -----------------------------
  !do while ( eps < 1d-1 )
  !  call RK4( q, z, dz, params, argsint, q_new )
  !  dz = dz*10.0_rk
  !  eps = maxval(abs(q(1,:)-q_new(1,:))/abs(q(1,:)))
  !end do
  !dz = dz / 10.0_rk

  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "starting main time loop"
  end if

  atol = 1e-10_rk ! absolute tolerance
  rtol = 1e-10_rk ! relative tolerance

  !options = set_opts( DENSE_J = .true.,USER_SUPPLIED_JACOBIAN= .false.,&
                    !RELERR=rtol,ABSERR=atol )
                    !METHOD_FLAG=25,LOWER_BANDWIDTH=nrhs-1, UPPER_BANDWIDTH=nrhs-1,RELERR=rtol,ABSERR=atol )
  !options = set_opts( METHOD_FLAG=25, BANDED_J = .true.,USER_SUPPLIED_JACOBIAN= .false.,&
                    !  RELERR=rtol,ABSERR=atol, LOWER_BANDWIDTH=nrhs-1, UPPER_BANDWIDTH=nrhs-1 )
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL,H0=dz,HMAX=0.0001_rk,MXSTEP=100000)
  itask = 1
  istate = 1
  write(*,*) params%gaxx
  write(*,*) params%gaff
  allocate(Y(nrhs*params%N))
  Y = reshape(q, (/nrhs*params%N/) )
  write(*,*) Y, z
  it = 1
  do while ( z <= params%z_max)!params%z_max .or. it < params%nt)
    !Y = pack(q, .true.)
    it = it + 1
    !zpdz=1.0_rk
    zpdz = z + params%dz_plot
    CALL VODE_F90(region3a_log,nrhs*params%N,Y,z,zpdz,itask,istate,OPTIONS,params,argsint)
    q_new = reshape(Y,(/nrhs, params%N/))
    z = zpdz
    !call RK4( q, z, dz, params, argsint, q_new )
    !q_tot(it,:,:) = q_new
    q_tot(it,1,:) = z
    q_tot(it,2:nrhs+1,:) = q_new
    T = params%mx/10**z
    q_tot(it,nrhs+2,:) = neq(T, params%mx, gDM)/ent(T, params)
    do i=1,params%N
      q_tot(it,nrhs+3,i) = neq(sqrt(params%gaff(i))*Ta(T,params), params%mx, gDM)/ent(T, params)
    end do
    !write(*,*) Y, z

  !  z = z + dz
    !write(*,*) "q=", q, z
  end do
  do i=1,params%N
    call write_matrix("temp/"//float2str(params%gaff(i))//float2str(params%gaxx(i))//".txt",q_tot(:,:,i))
  end do
!Y_fin = q_tot(end,1,:);
!do i=1,
!  if params%DM_low/params%mx < Y_fin(i) .and. params%DM_low/params%mx < Y_fin(i)
!    write(*,*) params%gaff(i), params%gaxx(i)
!  end if
!end do
  call MPI_Finalize(ierr)
end program main
