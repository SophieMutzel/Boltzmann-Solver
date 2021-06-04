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
  real(kind=rk)                        :: rtol, atol, helper, T, abserr,s, gam_xxff, sv_aaxx, H, sv_xxaa
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
  nrhs = 2
  call get_params( params )
  ! create parameter struct
  call allocate_couplings( params )

  call ini_cons_to_params( params, argsint )

  argsint%mx = params%mx
  argsint%ma = params%ma
  allocate(q(nrhs,params%N))
  L = int((params%z_max-params%z_start)/params%dz_plot)
  allocate(q_tot(2*nrhs+1,params%N,L+3))
  allocate(rhs(7,params%N,L+3))
  T = params%mx/10**params%z_start
  s = ent(T,params)
  do i=1,params%N
    ! Y_chi=Yeq(T')
    q(1,i) = neq(T, params%mx, gDM)
    q(2,i) = neq(T, params%ma, ga)
  end do
  q_tot(1,:,1) = params%z_start
  q_tot(2,:,1) = q(1,:)/s
  q_tot(3,:,1) = q(2,:)/s
  q_tot(4,:,1) = q(1,:)/s
  q_tot(5,:,1) = q(2,:)/s

!  ! check if equilibrium in HS attained
!  call sigmav( T, params, argsint, "aaxx", test )
!  if (test*params%gaxx(1)*params%gaxx(1)*params%gaxx(1)&
!      *params%gaxx(1)*neq(T, params%ma, ga)<Hub(params%mx/10**params%z_start,0.0_rk ))&
!      write(*,*) "Eq not attained!"

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
    call rhs_contributions_in_n( nrhs*params%N, z, Y, params, argsint, rhs(:,:,it), "fo" )
    it = it + 1
    zpdz = z + params%dz_plot

    CALL VODE_F90( region_freeze_out_cp, nrhs*params%N, Y, z, zpdz, itask, &
                   istate, OPTIONS, params, argsint )
    q = reshape(Y,(/nrhs, params%N/))
    z = zpdz
    T = params%mx/10**z
    s = ent(T,params)
    H = Hub(T,rhoeq(T,params%ma,ga)+rhoeq(T,params%mx,gDM))
    q_tot(1,:,it) = z
    q_tot(2:nrhs+1,:,it) = q/s
    q_tot(4,:,it) = neq(T, params%mx, gDM)/s
    q_tot(5,:,it) = neq(T, params%ma, ga)/s

    eps = abs((q_tot(1,1,it)-q_tot(1,1,it-1))/q_tot(1,1,it-1))

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
    call write_gnuplot(98, "temp/"//trim(adjustl(params%file)),Yxmxobs/params%mx, "fo")
    call write_matrix("temp/rhs_"//trim(adjustl(params%file))//".txt",rhs(:,i,1:it))
    call write_gnuplot_rhs(99, "temp/rhs_"//trim(adjustl(params%file)))
!    if (io_error==0) then
!      write(97,*) params%gaxx(i)*params%gaff(i), params%gaxx(i), q_tot(2,i,it)
!    else
!      write(*,*) 'error', io_error,' while opening the file temp/in_n.txt'
!    end if
  end do
  !close(97)

  call MPI_Finalize(ierr)
end program main
