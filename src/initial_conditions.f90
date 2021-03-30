subroutine initial_conditions( params, q, q_new, q_tot, rhs )
  implicit none
  type (type_params), intent(in)                                 :: params
  real(kind=rk), dimension(:,:), allocatable, intent(inout)      :: q, q_new
  real(kind=rk), dimension(:,:,:), allocatable, intent(inout)    :: q_tot, rhs
  integer(kind=ik)                                               :: L, i
  real(kind=rk)                                                  :: T_start, s

  allocate(q(nrhs,params%N), q_new(nrhs,params%N))
  L = int((params%z_max-params%z_start)/params%dz_plot)
  allocate(q_tot(nrhs+5,params%N,L+3))
  allocate(rhs(nrhs+5,params%N,L+3))
  !allocate(q_tot(nrhs+5,params%N,params%nt/100))
  T_start = params%mx/10**params%z_start
  s = ent(T_start, params)
  do i=1,params%N
    ! Y_chi=Yeq(T')
    q(1,i) = neq(Ta(T_start,params), params%mx, gDM)/s
    !q(1,i) = neq(sqrt(params%gaff(i))*Ta(T_start,params), params%mx, gDM)/s
    ! Y_a=Yeq,a(T')
    q(2,i) = neq(Ta(T_start,params), params%ma, ga)/s
    !q(2,i) = neq(sqrt(params%gaff(i))*Ta(T_start,params), params%ma, ga)/s
  end do
  ! Y_chi
  !q(1,:) = params%initial_values(1) !* params%gaxx*params%gaxx*params%gaxx*params%gaxx
  ! Y_a
  !q(2,:) = params%initial_values(2) !* params%gaff*params%gaff
  !rho_a/rho
  !q(3,:) = params%initial_values(3) * params%gaff*params%gaff

  q_new = q
  q_tot(1,:,1) = params%z_start
  q_tot(2:nrhs+1,:,1) = q
  q_tot(nrhs+2,:,1) = neq(T_start, params%mx, gDM)/s
  q_tot(nrhs+3,:,1) = q(1,:)
  q_tot(nrhs+4,:,1) = q(2,:)
  q_tot(nrhs+5,:,1) = params%mx/Ta(T_start,params)!params%mx/sqrt(params%gaff)/Ta(T_start,params)

end subroutine initial_conditions
