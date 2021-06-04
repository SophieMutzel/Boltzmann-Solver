! thermal gluon mass
real(kind=rk) function mg_th(T)
  implicit none
  real(kind=rk), intent(in)     :: T
  mg_th = sqrt(4.0_rk*pi*alpha_s(2.0_rk*pi*T)*T*T*(3.0_rk + dble(count(mf(4:9)<T))/2.0_rk)/3.0_rk)
  return
end function mg_th

! thermal photon mass
real(kind=rk) function mgamma_th(T)
  implicit none
  real(kind=rk), intent(in)     :: T

  mgamma_th = sqrt(4.0_rk*pi*alpha_qed_th(2.0_rk*pi*T)*T*T*dble(count(mf<T))/3.0_rk)
end function mgamma_th

! thermal quark mass
real(kind=rk) function mq_th(T,mfi,qf)
  implicit none
  real(kind=rk), intent(in)     :: T,mfi,qf
  real(kind=rk)                 :: mq
  if (T>QCDcut) then
    mq = 4.0_rk*pi*alpha_s(2.0_rk*pi*T)*T*T/6.0_rk
  else
    mq = 0
  end if
  mq_th = mfi + Sqrt(mq+4.0_rk*pi*alpha_qed_th(2.0_rk*pi*T)*T*T*qf*qf/8.0_rk)
  return
end function mq_th

! thermal lepton mass
real(kind=rk) function ml_th(T,mfi,qf)
  implicit none
  real(kind=rk), intent(in)     :: T,mfi,qf
  ml_th = mfi + Sqrt(4.0_rk*pi*alpha_qed_th(2.0_rk*pi*T)*T*T*qf*qf/8.0_rk)
  return
end function ml_th
