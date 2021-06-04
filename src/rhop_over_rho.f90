real(kind=rk) function rhop_over_rho(LT,argsint)
  implicit none
  real(kind=rk), intent(in)         :: LT
  type (type_argsint), intent(in)   :: argsint
  real(kind=rk)                     :: drhoa,T
  integer(kind=ik)                  :: nd

  !T=10**LT
  T=LT
  nd = size(argsint%drhoa,2)
  call interp_linear(nd, argsint%drhoa(1,:),argsint%drhoa(2,:),log10(T), drhoa)
  ! inverse decays ff->a
  drhoa = drhoa + drho_decay(T, argsint%ma, "ffath")
  rhop_over_rho = -drhoa/(rho_SM(T)*Hub( T, 0.0_rk )*T)
end function rhop_over_rho
