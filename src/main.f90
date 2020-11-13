program main
!------------------------------------------------------------------------------!
  use mpi
  use module_params
  use module_utils
  use module_cosmo
  use module_xsecs
  use module_rhs
!------------------------------------------------------------------------------!
  implicit none
  ! global parameters
  real(kind=rk)                        :: z, dz, eps
  type (type_params)                   :: params
  type (type_argsint)                  :: argsint
  real(kind=rk), allocatable           :: q(:,:), q_new(:,:)
  real(kind=rk), allocatable           :: q_tot(:,:,:)
  integer(kind=ik)                     :: it, ierr, rank, nprocs

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
  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "optimizing dz"
  end if

  eps = 1d-3
  do while ( eps < 1d-1 )
    call RK4( q, z, dz, params, argsint, q_new )
    dz = dz*10.0_rk
    eps = maxval(abs(q(1,:)-q_new(1,:))/abs(q(1,:)))
  end do
  dz = dz / 10.0_rk
  
  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "starting main time loop"
  end if

  do while ( z < params%z_max .or. it < params%nt)
    it = it + 1
    call RK4( q, z, dz, params, argsint, q_new )
    q_tot(it,:,:) = q_new
    q = q_new
    z = z + dz
  end do

!Y_fin = q_tot(end,1,:);
!do i=1,
!  if params%DM_low/params%mx < Y_fin(i) .and. params%DM_low/params%mx < Y_fin(i)
!    write(*,*) params%gaff(i), params%gaxx(i)
!  end if
!end do
  call MPI_Finalize(ierr)
end program main
