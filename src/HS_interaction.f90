real(kind=rk) function kernel_aaxx( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(in)     :: argsint
    kernel_aaxx = sigma_aaxx(s,argsint%mx,argsint%ma,argsint%g)*4.0_rk*&
                  F(s,argsint%ma,argsint%ma)*F(s,argsint%ma,argsint%ma)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
    return
end function kernel_aaxx

real(kind=rk) function kernel_xxaa( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(in)     :: argsint

    kernel_xxaa = sigma_xxaa(s,argsint%mx,argsint%ma,argsint%g)*4.0_rk*&
                  F(s,argsint%mx,argsint%mx)*F(s,argsint%mx,argsint%mx)/sqrt(s) * bessK1( sqrt(s)/argsint%T)

end function kernel_xxaa

real(kind=rk) function kernel_axax( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(in)     :: argsint

    kernel_axax = sigma_axax(s,argsint%mx,argsint%ma,argsint%g)*4.0_rk*&
                  F(s,argsint%mx,argsint%ma)*F(s,argsint%mx,argsint%ma)/sqrt(s) * bessK1( sqrt(s)/argsint%T)

end function kernel_axax

real(kind=rk) function F(s,m1,m2)
  implicit none
    real(kind=rk), intent(in)  :: s, m1, m2
    F = 0.5_rk*sqrt((s-(m1+m2)*(m1+m2))*(s-(m1-m2)*(m1-m2)))
    return
end function F

real(kind=rk) function kernel_xxaa_series( s, argsint )
  implicit none
  real(kind=rk), intent(in)           :: s
  type (type_argsint), intent(in)     :: argsint
  real(kind=rk)                       :: mx, T

  mx = argsint%mx
  T = argsint%T

  kernel_xxaa_series = sigma_xxaa(s,mx,argsint%ma,argsint%g)*&
                F(s,mx,mx)*F(s,mx,mx)/sqrt(s)*&
                (exp((2.0_rk*mx - sqrt(s))/T)* mx* &
                (1140.0_rk*s*T*T -60.0_rk*mx*T*(8.0_rk*s + 3.0_rk*sqrt(s)*T)&
                +mx*mx*(128.0_rk*s + 48.0_rk*sqrt(s)*T - 15.0_rk*T*T)))&
                /(128.0_rk*mx*mx*mx*mx*mx*sqrt(2.0_rk*pi)*s**1.25_rk*T**1.5_rk)

end function kernel_xxaa_series

real(kind=rk) function kernel_aaxx_series( s, argsint )
  implicit none
  real(kind=rk), intent(in)           :: s
  type (type_argsint), intent(in)     :: argsint
  real(kind=rk)                       :: ma, T

  ma = argsint%ma
  T = argsint%T

  kernel_aaxx_series = sigma_aaxx(s,argsint%mx,ma,argsint%g)*&
                F(s,ma,ma)*F(s,ma,ma)/sqrt(s)*&
                (exp((2.0_rk*ma - sqrt(s))/T)* ma* &
                (1140.0_rk*s*T*T -60.0_rk*ma*T*(8.0_rk*s + 3.0_rk*sqrt(s)*T)&
                +ma*ma*(128.0_rk*s + 48.0_rk*sqrt(s)*T - 15.0_rk*T*T)))&
                /(128.0_rk*ma*ma*ma*ma*ma*sqrt(2.0_rk*pi)*s**1.25_rk*T**1.5_rk)

end function kernel_aaxx_series

real(kind=rk) function kernel_axax_series( s, argsint )
  implicit none
  real(kind=rk), intent(in)           :: s
  type (type_argsint), intent(in)     :: argsint
  real(kind=rk)                       :: mx, ma, T

  mx = argsint%mx
  ma = argsint%ma
  T = argsint%T

  kernel_axax_series = sigma_axax(s,mx,argsint%ma,argsint%g)*&
                F(s,mx,ma)*F(s,mx,ma)/sqrt(s)*&
                exp((ma + mx - Sqrt(s))/T)*(1.0_rk/((ma*mx*T)**1.5_rk*Sqrt(2.0_rk*Pi)*s**0.25)&
                + (3.0_rk*ma*mx - 15.0_rk*(ma + mx)*Sqrt(s))/&
                (8.0_rk*(ma*mx)**2.5*Sqrt(2.0_rk*Pi)*s**0.75*Sqrt(T))&
                + (15.0_rk*(-(ma*ma*mx*mx) - 6.0_rk*ma*mx*(ma + mx)*Sqrt(s) + &
                (23.0_rk*ma*ma + 30.0_rk*ma*mx + 23.0_rk*mx*mx)*s)*Sqrt(T))/&
                (128.*(ma*mx)**3.5*Sqrt(2.0_rk*Pi)*s**1.25))

end function kernel_axax_series
