subroutine initial_conditions( params, q, q_new, q_tot )
  implicit none
  type (type_params), intent(in)                                 :: params
  real(kind=rk), dimension(:,:), allocatable, intent(inout)      :: q, q_new
  real(kind=rk), dimension(:,:,:), allocatable, intent(inout)    :: q_tot
  integer(kind=ik)                                               :: L

  allocate(q(2,params%N), q_new(2,params%N))
  L = max(int(params%z_max/params%dz),params%nt)
  allocate(q_tot(L,2,params%N))

  q(1,:) = params%init_Y * params%gaff*params%gaff * params%gaxx*params%gaxx
  q(2,:) = params%init_rhoprime * params%gaff*params%gaff * params%gaxx*params%gaxx

  q_tot(1,:,:) = q

end subroutine initial_conditions
