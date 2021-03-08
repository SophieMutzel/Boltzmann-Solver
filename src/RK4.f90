subroutine RK4( u0, t, dt, params, argsint, u1 )
  implicit none
  real(kind=rk), intent(in)                        :: t, dt
  type (type_params), intent(in)                   :: params
  type (type_argsint), intent(inout)               :: argsint
  real(kind=rk), dimension(:,:), intent(in)        :: u0
  real(kind=rk), dimension(:,:), intent(out)       :: u1
  real(kind=rk), dimension(nrhs,params%N)          :: k1, k2, k3, k4
  !Runge Kutta 4. order

  !call rhs_boltzmann( u0, params, argsint, t, k1 )
  !call rhs_boltzmann( u0 + dt/2.0_rk* k1, params, argsint, t + dt/2.0_rk, k2 )
  !call rhs_boltzmann( u0 + dt/2.0_rk* k2, params, argsint, t + dt/2.0_rk, k3 )
  !call rhs_boltzmann( u0 + dt * k3, params, argsint, t + dt, k4 )

!  call rhs_region3a( u0, params, argsint, t, k1 )
!  call rhs_region3a( u0 + dt/2.0_rk* k1, params, argsint, t + dt/2.0_rk, k2 )
!  call rhs_region3a( u0 + dt/2.0_rk* k2, params, argsint, t + dt/2.0_rk, k3 )
!  call rhs_region3a( u0 + dt * k3, params, argsint, t + dt, k4 )
  call region3a_eq( params%N*nrhs, t, u0, k1, params, argsint )
  call region3a_eq( params%N*nrhs, t + dt/2.0_rk, u0 + dt/2.0_rk* k1, k2, params, argsint )
  call region3a_eq( params%N*nrhs, t + dt/2.0_rk, u0 + dt/2.0_rk* k2, k3, params, argsint )
  call region3a_eq( params%N*nrhs, t + dt, u0 + dt * k3, k4, params, argsint )

  u1 = u0 + dt/6.0_rk * ( k1 + 2.0_rk*k2 + 2.0_rk*k3 + k4)

end subroutine RK4
