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
  rhop_over_rho = -drhoa*argsint%g*argsint%g/(rho_SM(T)*Hub( T, 0.0_rk )*T)
end function rhop_over_rho

real(kind=rk) function n_axion_seq(T,argsint)
  implicit none
  real(kind=rk), intent(in)           :: T
  type (type_argsint), intent(inout)  :: argsint
  real(kind=rk)                       :: gam_agff, gam_afgf, ffa, s

  call gamma_r_new( T, argsint, "agffth", gam_agff )
  gam_afgf = 0.0_rk
  ffa = gammav(T, argsint, "affth")*argsint%g*argsint%g
  ! For now ignore HS degrees of freedom...
  s = 2.0_rk * pi* pi /45.0_rk * geff_s(T) * T*T*T
  n_axion_seq = l10/s/Hub(T,argsint%ra_ini)*(gam_agff + gam_afgf + ffa*neq(T, argsint%ma, ga))
end function n_axion_seq
