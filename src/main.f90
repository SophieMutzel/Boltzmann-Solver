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
  real(kind=rk)                        :: rtol, atol, helper, T, test, abserr, Tprime, s, ztest
  real(kind=rk), allocatable           :: Y(:), work(:), points(:), x_eval(:),y_eval(:)
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

  call initial_conditions( params, q, q_new, q_tot, rhs )

  argsint%mx = params%mx
  argsint%ma = params%ma

  dz = params%dz
  z  = params%z_start
  it = 0
  do while ( z <= params%z_max)
!    !call test_bezier(params%geff_HS(1,:),params%geff_HS(2,:),z,test,nd-1,params%A,params%B)
    zpdz = z +params%dz_plot
    T = params%mx/10**z
!    test = geff_rho(T)
!    call geffSM(T,params,Tprime)
    !Tprime = Ta(T,params)
    ! HS interaction
    !argsint%g = params%gaxx(1)
    !call sigmav( Tprime, params, argsint, "aaxx", test )
!    !test = 2.46743 - 0.900703 *z - 0.426853 *z**2 + 0.344933 *z**3 + 0.241269 *z**4 - 1.64352 *tanh(2.49447*z)
!    !call interp_linear(size(x_eval), x_eval,y_eval,z,test)!log10(Tprime), test)
    z = zpdz
    !write(*,*) z, test, Hub( T, params ), ent( T, params )
    call sigmav( T, params, argsint, "aaxx", test )
    write(*,*) z, T,  test
  end do
stop
  if (rank==0) then
    write(*,'(80("_"))')
    write(*,*) "starting main time loop"
  end if

  atol = 1e-13_rk ! absolute tolerance
  rtol = 1e-6_rk ! relative tolerance

  !options = set_opts( DENSE_J = .true.,USER_SUPPLIED_JACOBIAN= .false.,&
                    !RELERR=rtol,ABSERR=atol )
                    !METHOD_FLAG=25,LOWER_BANDWIDTH=nrhs-1, UPPER_BANDWIDTH=nrhs-1,RELERR=rtol,ABSERR=atol )
  !options = set_opts( METHOD_FLAG=25, BANDED_J = .true.,USER_SUPPLIED_JACOBIAN= .false.,&
                    !  RELERR=rtol,ABSERR=atol, LOWER_BANDWIDTH=nrhs-1, UPPER_BANDWIDTH=nrhs-1 )
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL,H0=dz,HMAX=0.0001_rk,MXSTEP=10000)
  itask = 1
  istate = 1
  allocate(Y(nrhs*params%N))
  Y = reshape(q, (/nrhs*params%N/) )
  it = 1

  call rhs_contributions_in_n( nrhs*params%N, z, Y, params, argsint, rhs(:,:,it) )

  eps = 1.0_rk
  conv_eps = 1e-3_rk
!  allocate(x_eval((nd-1)*4),y_eval((nd-1)*4))
!  call test_bezier(params%geff_HS(1,:),params%geff_HS(2,:),x_eval,y_eval,nd-1,4,params%A,params%B)
!  do i=1,(nd-1)*4
!    write(*,*) x_eval(i), y_eval(i)
!  end do
!stop
  z  = params%z_start
  do while ( z <= params%z_max .and. eps > conv_eps .or. z<=0.0_rk)
    it = it + 1
    zpdz = z + params%dz_plot
    call rhs_contributions_in_n( nrhs*params%N, z, Y, params, argsint, rhs(:,:,it) )

    CALL VODE_F90( region3a_in_n, nrhs*params%N, Y, z, zpdz, itask, &
                   istate, OPTIONS, params, argsint )
    q_new = reshape(Y,(/nrhs, params%N/))
    z = zpdz
    T = params%mx/10**z
    s = ent(T,params)
    q_tot(1,:,it) = z
    q_tot(2:nrhs,:,it) = q_new(1:2,:)/s
    !Tprime=Tanew(T,params,q_new(3,:),q_new(1,:)/neq(T,params%mx,gDM)*s*rhoeq(T,params%mx,gDM))
    !Tprime = Ta(T,params)
    Tprime = q_new(3,1)
    q_tot(nrhs+1,:,it) = neq(T, params%mx, gDM)/s
    do i=1,params%N
      q_tot(nrhs+2,i,it) = neq(Tprime, params%mx, gDM)/s
      q_tot(nrhs+3,i,it) = neq(Tprime, params%ma, ga)/s
    end do
    q_tot(nrhs+4,:,it) = Tprime!params%mx/(sqrt(params%gaff)*Tprime)

    eps = max(abs((q_tot(2,1,it)-q_tot(2,1,it-1))/q_tot(2,1,it-1)),abs((q_tot(3,1,it)-q_tot(3,1,it-1))/q_tot(3,1,it-1)))
    !if (z>1.0_rk ) write(*,*) z, q_new(1,1)/(exp(-params%mx/Tprime)*(gDM*(params%mx*Tprime/pi)**(1.5_rk)/(2.0_rk* sqrt(2.0_rk)))),&
    !                q_new(2,1)/(exp(-params%ma/Tprime)*(ga*(params%ma*Tprime/pi)**(1.5_rk)/(2.0_rk* sqrt(2.0_rk)))), Tprime
    write(*,*) z, eps, Tprime
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
  open (unit=97, file="temp/in_n.txt", status='old', action='write', position='append', iostat=io_error)
  do i=1,params%N
  !  call write_matrix("temp/"//exp2str(params%gaff(i))//exp2str(params%gaxx(i))//".txt",q_tot(:,i,:))
    call write_matrix("temp/"//trim(adjustl(params%file))//".txt",q_tot(:,i,1:it))
    call write_matrix("temp/rhs_"//trim(adjustl(params%file))//".txt",rhs(:,i,1:it))
    call write_gnuplot(98, "temp/"//trim(adjustl(params%file)),Yxmxobs/params%mx)
    call write_gnuplot_rhs(99, "temp/rhs_"//trim(adjustl(params%file)))
    if (io_error==0) then
      write(97,*) params%gaxx(i)*params%gaff(i), params%gaxx(i), q_tot(2,i,it)
    else
      write(*,*) 'error', io_error,' while opening the file temp/in_n.txt'
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
