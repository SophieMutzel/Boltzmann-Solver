module module_xsecs

use module_precision
use mpi

contains
  real(kind=rk) function sigma_xxaa(s, mx, ma )
    implicit none
    real(kind=rk), intent(in)        :: s, mx, ma
    real(kind=rk)                    :: mxsq, masq

    masq = ma*ma
    mxsq = mx*mx
    sigma_xxaa = 0.5_rk*( mxsq*mxsq*((Sqrt(-4.0_rk*masq+s)*&
                (-3.0_rk*masq*masq+8.0_rk*masq*mxsq-2.0_rk*mxsq*s))&
                /(Sqrt(-4.0_rk*mxsq+s)* (masq*masq-4.0_rk*masq*mxsq+mxsq*s))+((6.0_rk*masq*masq-4.0_rk*masq*s+s*s)*&
                (Log(-2.0_rk*masq+s-Sqrt((-4.0_rk*masq+s)*(-4.0_rk*mxsq+s)))-Log(-2.0_rk*masq+s+Sqrt((-4.0_rk*masq+s)*&
                (-4.0_rk*mxsq+s)))))/((2.0_rk*masq-s)*(-4.0_rk*mxsq+s))))/(8.0_rk*pi*s)


    ! 1/(8 \[Pi] s (-4 ma^2+s)) gaxx^4 mx^4 (-4 mx^2+s) ((Sqrt[-4 ma^2+s] (-3 ma^4+8 ma^2 mx^2-2 mx^2 s))/
    !(Sqrt[-4 mx^2+s] (ma^4-4 ma^2 mx^2+mx^2 s))+1/((2 ma^2-s) (-4 mx^2+s)) (6 ma^4-4 ma^2 s+s^2)
    !(Log[-2 ma^2+s-Sqrt[(-4 ma^2+s) (-4 mx^2+s)]]-Log[-2 ma^2+s+Sqrt[(-4 ma^2+s) (-4 mx^2+s)]]))
    return
  end function sigma_xxaa

  real(kind=rk) function sigma_aaxx(s, mx, ma )
    implicit none
    real(kind=rk), intent(in)        :: s, mx, ma
    real(kind=rk)                    :: mxsq, masq

    masq = ma*ma
    mxsq = mx*mx

    sigma_aaxx = 4.0_rk*(mxsq*mxsq*((Sqrt(-4.0_rk*mxsq+s)*(-3.0_rk* masq*masq+8.0_rk*masq*mxsq-2.0_rk* mxsq*s))&
                /(Sqrt(-4.0_rk* masq+s)*(masq*masq-4.0_rk* masq*mxsq+mxsq*s))&
                -((6.0_rk* masq*masq-4.0_rk* masq*s+s*s)*&
                (Log(-2.0_rk* masq+s-Sqrt((-4.0_rk* masq+s)*(-4.0_rk* mxsq+s)))&
                -Log(-2.0_rk* masq+s+Sqrt((-4.0_rk* masq+s)&
                *(-4.0_rk* mxsq+s)))))/(8.0_rk* masq*masq-6.0_rk* masq*s+s*s)))/(8.0_rk*pi*s)
!    sigma_aaxx = 1.0_rk/(8.0_rk*pi*s*(-4.0_rk*masq+s))*mxsq*mxsq*(-4.0_rk*mxsq+s)*&
!                ((sqrt(-4.0_rk*masq+s)*(-3.0_rk*masq*masq+8.0_rk*masq*mxsq-2.0_rk*mxsq*s))/(Sqrt(-4.0_rk*mxsq+s)*&
!                (masq*masq-4.0_rk*masq*mxsq+mxsq*s))+1.0_rk/((2.0_rk*masq-s) *&
!                (-4.0_rk*mxsq+s))*(6.0_rk*masq*masq-4.0_rk*masq*s+s*s)*&
!                (Log(-2.0_rk*masq+s-Sqrt((-4.0_rk*masq+s)*(-4.0_rk*mxsq+s)))&
!                -Log(-2.0_rk*masq+s+Sqrt((-4.0_rk*masq+s)*(-4.0_rk*mxsq+s)))))
    return
  end function sigma_aaxx

  real(kind=rk) function sigma_xxff(s, mf, nc, mx, ma)
    implicit none
    real(kind=rk), intent(in)        :: s, mf, nc, mx, ma


    sigma_xxff = nc * s * sqrt(-4.*mf*mf + s) * mf*mf / &
                ( 16. * pi * (ma*ma - s)*(ma*ma - s) * real( sqrt(-4.*mx*mx + s) ) )
    return
  end function sigma_xxff

  real(kind=rk) function sigma_ffa(s, mf, ma)
    implicit none
    real(kind=rk), intent(in)         :: s, mf, ma
    real(kind=rk)                     :: Gamma

    Gamma = ma * mf*mf /8.0_rk/pi * sqrt( 1.0_rk - 4.0_rk*mf*mf/ma/ma )

    sigma_ffa = Gamma/( (s-ma*ma)*(s-ma*ma) + ma*ma*Gamma*Gamma) *sqrt(s)*mf*mf/2.0_rk
    return
  end function sigma_ffa

  ! prefac for gamma = alpha
  ! prefac for g = alpha_s*4/9
  real(kind=rk) function sigma_agff(s, mf, ma)
    implicit none
    real(kind=rk), intent(in)         :: s, mf, ma
    real(kind=rk)                     :: masq, mfsq

    masq = ma*ma
    mfsq = mf*mf
    !sigma_ffga = prefac* ( mfsq*qf*qf* (masq*sqrt((4.0_rk mfsq-s)*s) - (masq*masq - 4.0_rk*masq*mfsq + s*s)&
    !                atan(sqrt(-1.0_rk+(4.0_rk*mfsq)/s))))/(2.0_rk*pi*sqrt((masq-s)*(masq-s))*sqrt(-s*s*(-4.0_rk mfsq+s)*(-4.0_rk*mfsq+s)))
    sigma_agff = 2.0_rk*( mfsq* sqrt(s* (-4.0_rk*mfsq+s))*(-2.0_rk*masq*sqrt(1.0_rk-(4.0_rk* mfsq)/s)*s+(masq*masq-4.0_rk*masq*mfsq+s*s)* &
                  log((-2.0_rk*mfsq+s+sqrt(1.0_rk-(4.0_rk*mfsq)/s)*s)/(2.0_rk*mfsq))))/(4.0_rk*pi*sqrt(1.0_rk-(4.0_rk*mfsq)/s)*((masq-s)*(masq-s))**(1.5)*s)
    return
  end function sigma_agff

  real(kind=rk) function sigma_afgf(s, mf, ma)
    implicit none
    real(kind=rk), intent(in)         :: s, mf, ma
    real(kind=rk)                     :: masq, mfsq

    masq = ma*ma
    mfsq = mf*mf
    sigma_afgf = 2.0_rk*(mfsq*(mfsq*mfsq-2.0_rk*mfsq*s+s*s)*sqrt(masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s))* &
                ((-mfsq+s)*(-mfsq+s)*(-mfsq+3.0_rk*s)+masq*(mfsq*mfsq-2.0_rk*mfsq*s-7.0_rk*s*s)-(2.0_rk*s*s* &
                (2.0_rk*masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s))* &
                (-log(-masq+mfsq+s-sqrt(masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s)))&
                +log(-masq+mfsq+s+sqrt(masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s)))))&
                /sqrt(masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s))))/&
                (16.0_rk*pi*(mfsq-s)*(mfsq-s)*(mfsq-s)*s*s*&
                (masq*masq-2.0_rk*masq*mfsq+mfsq*mfsq-2.0_rk*masq*s-2.0_rk*mfsq*s+s*s))
    return
end function sigma_afgf

end module
