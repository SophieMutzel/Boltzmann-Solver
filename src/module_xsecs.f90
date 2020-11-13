module module_xsecs

use module_precision
use mpi

contains
  real(kind=rk) function sigma_xxaa(s, mf, nc, mx, ma )
    implicit none
    real(kind=rk), intent(in)        :: s, mf, nc, mx, ma

    sigma_xxaa = -1./( 8. * pi * s * sqrt(-4.*mx*mx + s) )&
                * ( -4. *ma*ma + s)**2 * sqrt(s - 4.*ma*ma) *(-4.*mx*mx + s) &
                * ( (3.*ma**4 - 8.*ma*ma*mx*mx + 2*mx*mx*s)/&
                ( (-4. *ma*ma + s)**2 *(-4.*mx*mx + s) *(ma**4 - 4.*ma*ma*mx*mx + mx*mx*s)) &
                -( (6.*ma**4 - 4.*ma*ma*s + s*s) * (log(-2.*ma*ma + s - sqrt((-4.*ma*ma + s)*(-4.*mx*mx + s))) &
                -log(-2.*ma*ma + s + sqrt((-4.*ma*ma + s)*(-4.*mx*mx + s)))))/ &
                ((2.*ma*ma - s)*(-4.*ma*ma + s)**(2.5) *(-4.*mx*mx + s)**(1.5)))
    return
  end function sigma_xxaa

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
  end function sigma_ffa

  ! prefac for gamma = alpha
  ! prefac for g = alpha_s*4/9
  real(kind=rk) function sigma_ffga(s, mf, ma, qf, prefac)
    implicit none
    real(kind=rk), intent(in)         :: s, mf, ma, qf, prefac
    real(kind=rk)                     :: masq, mfsq

    masq = ma*ma
    mfsq = mf*mf
    sigma_ffga = prefac* ( mfsq*qf*qf* (masq*sqrt((4.0_rk mfsq-s)*s) - (masq*masq - 4.0_rk*masq*mfsq + s*s)&
                    atan(sqrt(-1.0_rk+(4.0_rk*mfsq)/s))))/(2.0_rk*pi*sqrt((masq-s)*(masq-s))*sqrt(-s*s*(-4.0_rk mfsq+s)*(-4.0_rk*mfsq+s)))
  end function sigma_ffga

  ! prefac for gamma = alpha
  !prefac for g = alpha_s/6
  real(kind=rk) function sigma_fgfa(s, mf, ma, qf, prefac)
    implicit none
    real(kind=rk), intent(in)         :: s, mf, ma, qf, prefac
    real(kind=rk)                     :: masq, mfsq

    masq = ma*ma
    mfsq = mf*mf
    sigma_fgfa = prefac*(1.0_rk/(64.0_rk*pi*((mfsq-s)*(mfsq-s)*(mfsq-s))*s*s))*mfsq*qf*qf*sqrt(masq*masq+(mfsq-s)*(mfsq-s)-2.0_rk*masq*(mfsq+s))*&
                (-(mfsq-3.0_rk*s)*(mfsq-s)*(mfsq-s)+masq*(mfsq*mfsq-2.0_rk*mfsq*s-7.0_rk*s*s)+2.0_rk*s*s*(-(mfsq-s)*(mfsq-s)+2.0_rk*masq*(mfsq+s)+&
                masq*masq*(-1.0_rk-1.0_rk/(masq*masq+(mfsq-s)*(mfsq-s)-2.0_rk*masq*(mfsq+s))))*log(-1+(2.0_rk*(masq-mfsq-s))&
                  /(masq*masq+mfsq*mfsq+masq*(1.0_rk-2.0_rk*mfsq-2.0_rk*s)+(-1.0_rk+s)*s-mfsq*(1.0_rk+2.0_rk*s))))
  end function sigma_fgfa

end module
