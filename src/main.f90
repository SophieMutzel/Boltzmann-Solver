program main
  use module_params
  use module_utils
  use module_cosmo
  use module_xsecs
  use module_rhs
  implicit none
  ! global parameters
  real(kind=rk)                        :: z, dz, gHS, result
  type (type_params)                   :: params
  type (type_argsint)                  :: argsint
  real(kind=rk), allocatable           :: q(:,:), q_new(:,:)
  real(kind=rk), allocatable           :: q_tot(:,:,:)
  integer(kind=ik)                     :: it

  ! create parameter struct
  call allocate_couplings( params )

  call ini_cons_to_params( params )

  call initial_conditions( params, q, q_new, q_tot )

  argsint%mx = params%mx
  argsint%ma = params%ma
  dz = params%dz
  z = params%z_start
  it = 0
  write(*,*) "starting main time loop"
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

end program main
