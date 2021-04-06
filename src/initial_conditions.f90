subroutine initial_conditions( params, q, q_new, q_tot, rhs )
  implicit none
  type (type_params), intent(in)                                 :: params
  real(kind=rk), dimension(:,:), allocatable, intent(inout)      :: q, q_new
  real(kind=rk), dimension(:,:,:), allocatable, intent(inout)    :: q_tot, rhs
  integer(kind=ik)                                               :: L, i, nr
  real(kind=rk)                                                  :: T_start, s, rar, rho,Tprime

  allocate(q(nrhs,params%N), q_new(nrhs,params%N))
  L = int((params%z_max-params%z_start)/params%dz_plot)
  allocate(q_tot(nrhs+5,params%N,L+3))
  allocate(rhs(nrhs+5,params%N,L+3))
  !allocate(q_tot(nrhs+5,params%N,params%nt/100))
  T_start = params%mx/10**params%z_start
!  nr = size(params%rhoa_rho,2)
!  call interp_linear(nr, params%rhoa_rho(1,:),params%rhoa_rho(2,:),T_start, rar)
!  rho = rar*geff_rho(T_start)*pi*pi/30.0_rk*T_start*T_start*T_start*T_start
  s = ent(T_start, params)
  !Tprime = Tanew(T_start,params,(/rho/),(/1.0_rk/))
  Tprime = Ta(T_start,params)
  do i=1,params%N
    ! Y_chi=Yeq(T')
    q(1,i) = neq(Tprime, params%mx, gDM)/s
    ! Y_a=Yeq,a(T')
    q(2,i) = neq(Tprime, params%ma, ga)/s
    ! rho'=rhoeq,x(T')+rhoeq,a(T')
    q(3,i) = Tprime!rhoeq(Tprime, params%ma, ga)+rhoeq(Tprime, params%mx, gDM)
  end do

  q_new = q
  q_tot(1,:,1) = params%z_start
  q_tot(2:nrhs,:,1) = q(1:2,:)
  q_tot(nrhs+1,:,1) = neq(T_start, params%mx, gDM)/s
  q_tot(nrhs+2,:,1) = q(1,:)
  q_tot(nrhs+3,:,1) = q(2,:)
  q_tot(nrhs+4,:,1) = Tprime

end subroutine initial_conditions
