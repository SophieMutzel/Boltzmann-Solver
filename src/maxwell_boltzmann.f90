real(kind=rk) function neq( T, m, g )
  implicit none
  real(kind=rk), intent(in)           :: T, m, g
    neq = 0.5_rk * g * T * m*m * bessK2(m/T) /pi/pi
end function neq

real(kind=rk) function rhoeq( T, m, g )
  implicit none
  real(kind=rk), intent(in)           :: T, m, g

  rhoeq = 0.5_rk * g *m*m*T*(m*bessK1(m/T)+3.0_rk*T*bessK2(m/T)) /pi/pi
end function rhoeq

real(kind=rk) function peq( T, m, g )
  implicit none
  real(kind=rk), intent(in)           :: T, m, g

  peq = 0.5_rk * g *m*m*T*T*bessK2(m/T)/pi/pi

end function peq

real(kind=rk) function drhoeq( T, m, g )
  ! d(rho_eq(T)/dT
  implicit none
  real(kind=rk), intent(in)           :: T, m, g
  drhoeq = g*m/(2.0_rk*pi*pi)*((m*m*m/T + 12.0_rk*m*T)*bessK0(m/T) + (5.0_rk*m*m + 24.0_rk*T*T)*bessK1(m/T))
end function drhoeq

real(kind=rk) function rhoeqneq( T, m )
  ! rho_eq(T)/neq(T)
  implicit none
  real(kind=rk), intent(in)           :: T, m
  if (T/m>0.03_rk) then
    rhoeqneq = 3*T + m*bessK1(m/T)/bessK2(m/T)
  else
    rhoeqneq = m + (3.0_rk*T)/2.0_rk+ (15.0_rk*T*T)/(8.0_rk*m)
  end if
end function rhoeqneq
real(kind=rk) function drhoeqneq( T, m )
  ! d(rho_eq(T)/neq(T))/dT
  implicit none
  real(kind=rk), intent(in)           :: T, m
  real(kind=rk)                       :: bk1, bk2, T2,m2
  if (T/m<0.03_rk) then
    T2 = T*T
    m2 = m*m
    drhoeqneq = 3.0_rk/2.0_rk + (15.0_rk*T)/(4.0_rk*m) - &
                (45.0_rk*T2)/(8.0_rk*m2) +((135.0_rk*T*T2)/(32.0_rk*m*m2))&
                +((225.0_rk*T2*T2)/(32.0_rk*m2*m2))-&
                ((22275.0_rk*T2*T2*T)/(512.0_rk*m2*m2*m))+&
                (4725.0_rk*T2*T2*T2)/(32.0_rk*m2*m2*m2)&
                -((1905525.0_rk*T2*T2*T2*T)/(4096.0_rk*m2*m2*m2*m))
  else
  bk1 = bessK1(m/T)
  bk2 = bessK2(m/T)
  drhoeqneq = -((m*m* bk1*bk1 - bk2*(m*m*bessK0(m/T) + (m*m + 6.0_rk* T*T)*bk2) + &
              m*m*bk1*bessk_s(3,m/T))/(2.0_rk* T*T* bk2*bk2))
  end if
end function drhoeqneq
