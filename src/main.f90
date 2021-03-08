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
  real(kind=rk), allocatable           :: q(:,:), q_new(:,:)
  real(kind=rk), allocatable           :: q_tot(:,:,:), rhs(:,:,:)
  integer(kind=ik)                     :: it, ierr, rank, nprocs, N,itask
  integer(kind=ik)                     :: io_error, istate, i, ier, neval, nit, nd
  real(kind=rk)                        :: rtol, atol, helper, T, test, abserr
  real(kind=rk), allocatable           :: Y(:), work(:), points(:)
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

  call get_params( params )
  ! create parameter struct
  call allocate_couplings( params )

  call ini_cons_to_params( params, argsint )

  params%mx=10.0_rk
  call initial_conditions( params, q, q_new, q_tot, rhs )

  argsint%mx = params%mx
  argsint%ma = params%ma

  dz = params%dz
  z  = params%z_start
  it = 0

  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "starting main time loop"
  end if

  atol = 1e-18_rk ! absolute tolerance
  rtol = 1e-12_rk ! relative tolerance

  !options = set_opts( DENSE_J = .true.,USER_SUPPLIED_JACOBIAN= .false.,&
                    !RELERR=rtol,ABSERR=atol )
                    !METHOD_FLAG=25,LOWER_BANDWIDTH=nrhs-1, UPPER_BANDWIDTH=nrhs-1,RELERR=rtol,ABSERR=atol )
  !options = set_opts( METHOD_FLAG=25, BANDED_J = .true.,USER_SUPPLIED_JACOBIAN= .false.,&
                    !  RELERR=rtol,ABSERR=atol, LOWER_BANDWIDTH=nrhs-1, UPPER_BANDWIDTH=nrhs-1 )
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL,H0=dz,HMAX=0.0001_rk,MXSTEP=100000)
  itask = 1
  istate = 1
  do it = 1, size(params%geff_HS(2,:))
    write(*,*) params%geff_HS(1,it), params%geff_HS(2,it)
  end do
  stop
  allocate(Y(nrhs*params%N))
  Y = reshape(q, (/nrhs*params%N/) )
  it = 1
  call rhs_contributions( nrhs*params%N, z, Y, params, argsint, rhs(:,:,it) )

  eps = 1.0_rk
  conv_eps = 1e-2_rk
  do while ( z <= params%z_max .and. eps > conv_eps)!params%z_max .or. it < params%nt)
    it = it + 1
    zpdz = z + params%dz_plot
    CALL VODE_F90( region3a_eq, nrhs*params%N, Y, z, zpdz, itask, &
                   istate, OPTIONS, params, argsint )
    q_new = reshape(Y,(/nrhs, params%N/))
    z = zpdz
    call rhs_contributions( nrhs*params%N, z, Y, params, argsint, rhs(:,:,it) )

    q_tot(1,:,it) = z
    q_tot(2:nrhs+1,:,it) = q_new
    T = params%mx/10**z
    q_tot(nrhs+2,:,it) = neq(T, params%mx, gDM)/ent(T, params)
    do i=1,params%N
      q_tot(nrhs+3,i,it) = neq(sqrt(params%gaff(i))*Ta(T,params), params%mx, gDM)/ent(T, params)
      q_tot(nrhs+4,i,it) = neq(sqrt(params%gaff(i))*Ta(T,params), params%ma, ga)/ent(T, params)
    end do
    q_tot(nrhs+5,:,it) = params%mx/(sqrt(params%gaff)*Ta(T,params))

    eps = max(abs((q_tot(2,1,it)-q_tot(2,1,it-1))/q_tot(2,1,it-1)),abs((q_tot(3,1,it)-q_tot(3,1,it-1))/q_tot(3,1,it-1)))
    write(*,*) z, eps

  end do
!  nit = 1
!  do while ( it < params%nt)
!    !Y = pack(q, .true.)
!    it = it + 1
!    !zpdz=1.0_rk
!    zpdz = z + params%dz_plot
!    !z = zpdz
!    !write(*,*) z
!    Y = reshape(q, (/nrhs*params%N/) )
!    if (mod(it,100)==0) then
!      call rhs_contributions( nrhs*params%N, z, Y, params, argsint )
!    end if
!    call RK4( q, z, dz, params, argsint, q_new )
!    !q_tot(it,:,:) = q_new
!    z = z + dz
!    if (mod(it,100)==0) then
!      q_tot(1,:,nit) = z
!      q_tot(2:nrhs+1,:,nit) = q_new
!      T = params%mx/10**z
!      q_tot(nrhs+2,:,nit) = neq(T, params%mx, gDM)/ent(T, params)
!      do i=1,params%N
!        q_tot(nrhs+3,i,nit) = neq(sqrt(params%gaff(i))*Ta(T,params), params%mx, gDM)/ent(T, params)
!        q_tot(nrhs+4,i,nit) = neq(sqrt(params%gaff(i))*Ta(T,params), params%ma, ga)/ent(T, params)
!        !q_tot(nrhs+4,i,it) = neq(T, params%ma, ga)/ent(T, params)
!      end do
!      q_tot(nrhs+5,:,nit) = params%mx/(sqrt(params%gaff)*Ta(T,params))
!      nit= nit+1
!    end if
!    q = q_new
!    !write(*,*) Y, z
!
!    !write(*,*) "q=", q, z
!  end do
  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "main time loop over, writing results to files"
  end if

  ! write results to file
  open (unit=97, file="temp/all.txt", status='old', action='write', position='append', iostat=io_error)
  do i=1,params%N
    !call write_matrix("temp/"//exp2str(params%gaff(i))//exp2str(params%gaxx(i))//".txt",q_tot(:,i,:))
    call write_matrix("temp/"//trim(adjustl(params%file))//".txt",q_tot(:,i,1:it))
    call write_matrix("temp/rhs_"//trim(adjustl(params%file))//".txt",rhs(:,i,1:it))
    call write_gnuplot(98, "temp/"//trim(adjustl(params%file)))
    call write_gnuplot_rhs(99, "temp/rhs_"//trim(adjustl(params%file)))
    if (io_error==0) then
      write(97,*) params%gaxx(i)*params%gaff(i), params%gaxx(i), q_tot(2,i,it)
    else
      write(*,*) 'error', io_error,' while opening the file temp/all.txt'
    end if
  end do
  close(97)


!Y_fin = q_tot(end,1,:);
!do i=1,
!  if params%DM_low/params%mx < Y_fin(i) .and. params%DM_low/params%mx < Y_fin(i)
!    write(*,*) params%gaff(i), params%gaxx(i)
!  end if
!end do
  call MPI_Finalize(ierr)
end program main
