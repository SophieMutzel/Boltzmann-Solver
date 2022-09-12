module module_xsecs

use module_precision
use mpi
use module_utils

contains
  real(kind=rk) function sigma_xxaa(s, mx, ma, gaxx )
    implicit none
    real(kind=rk), intent(in)        :: s, mx, ma, gaxx
    real(kind=rk)                    :: mxsq, masq

    masq = ma*ma
    mxsq = mx*mx

    sigma_xxaa = gaxx*gaxx*gaxx*gaxx*( mxsq*mxsq*((Sqrt(-4.0_rk*masq+s)*&
                (-3.0_rk*masq*masq+8.0_rk*masq*mxsq-2.0_rk*mxsq*s))&
                /(Sqrt(-4.0_rk*mxsq+s)* (masq*masq-4.0_rk*masq*mxsq+mxsq*s))+((6.0_rk*masq*masq-4.0_rk*masq*s+s*s)*&
                (Log(-2.0_rk*masq+s-Sqrt((-4.0_rk*masq+s)*(-4.0_rk*mxsq+s)))-Log(-2.0_rk*masq+s+Sqrt((-4.0_rk*masq+s)*&
                (-4.0_rk*mxsq+s)))))/((2.0_rk*masq-s)*(-4.0_rk*mxsq+s))))/(16.0_rk*pi*s)

    return
  end function sigma_xxaa

  real(kind=rk) function sigma_aaxx(s, mx, ma, gaxx )
    implicit none
    real(kind=rk), intent(in)        :: s, mx, ma, gaxx
    real(kind=rk)                    :: mxsq, masq

    masq = ma*ma
    mxsq = mx*mx

    sigma_aaxx = (gaxx*gaxx*gaxx*gaxx*mxsq*mxsq*(&
                (Sqrt(-4.0_rk*mxsq + s)*(-3.0_rk*masq*masq + 8.0_rk*masq*mxsq - 2.0_rk*mxsq*s))&
                /(Sqrt(-4.0_rk*masq + s)*(masq*masq - 4.0_rk*masq*mxsq + mxsq*s)) - &
                ((6.0_rk*masq*masq - 4.0_rk*masq*s + s*s)*(Log(-2.0_rk*masq + s&
                 - Sqrt((-4.0_rk*masq + s)*(-4.0_rk*mxsq + s))) - &
                 Log(-2.0_rk*masq + s + Sqrt((-4.0_rk*masq + s)*(-4.0_rk*mxsq + s)))))&
                 /(8.0_rk*masq*masq - 6.0_rk*masq*s + s*s)))/(2.0_rk*pi*s)

    return
  end function sigma_aaxx

  real(kind=rk) function sigma_axax(s, mx, ma, gaxx )
    implicit none
    real(kind=rk), intent(in)        :: s, mx, ma, gaxx
    real(kind=rk)                    :: mxsq, masq

    masq = ma*ma
    mxsq = mx*mx
    sigma_axax = gaxx*gaxx*gaxx*gaxx*mxsq*mxsq*2.0_rk*pi*&
                (((ma - mx)*(ma - mx) - s)*((ma + mx)*(ma + mx) - s)*&
                (-2*masq*mxsq*(mxsq - s)**4 + 2.*ma**10*(mxsq + s)&
                 + 4*ma**6*(mxsq - s)**2*(2.*mxsq + 3.*s) + &
                 mxsq*(mxsq - s)**4*(mxsq + 5.*s) +&
                  ma**8*(-7.*mxsq*mxsq + 4.*mxsq*s - 5.*s*s) - &
                  masq*masq*(mxsq - s)**2*(2.*mxsq*mxsq + 15.*mxsq*s + 3.*s*s)) &
                  - 4.*(2.*masq*masq + (mxsq - s)**2)*s*s*(-mxsq + s)*&
                  (-2.*masq - mxsq + s)*(-(masq - mxsq)**2 + mxsq*s)*&
                  atanh((masq*masq + (mxsq - s)**2 - 2.*masq*(mxsq + s))/&
                  ((masq - mxsq)**2 + 2.*masq*s - s*s)))/(64.*Pi*Pi*(mxsq - s)**2&
                  *(2*masq + mxsq - s)*s*s*((masq - mxsq)**2 - mxsq*s)&
                  *(masq*masq + (mxsq - s)**2 - 2.*masq*(mxsq + s)))
    return
  end function sigma_axax

  real(kind=rk) function sigma_xxff(s, mf, nc, mx, ma, gaffgaxx)
    implicit none
    real(kind=rk), intent(in)        :: s, mf, nc, mx, ma, gaffgaxx


    sigma_xxff = gaffgaxx*gaffgaxx*mf*mf*nc * s * sqrt(-4.*mf*mf + s) * mx * mx / &
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
  real(kind=rk) function sigma_agff(s, mf,nc, ma, gaff)
    implicit none
    real(kind=rk), intent(in)         :: s, mf, ma, gaff, nc
    real(kind=rk)                     :: masq, mfsq

    masq = ma*ma
    mfsq = mf*mf
    !sigma_ffga = prefac* ( mfsq*qf*qf* (masq*sqrt((4.0_rk mfsq-s)*s) - (masq*masq - 4.0_rk*masq*mfsq + s*s)&
    !                atan(sqrt(-1.0_rk+(4.0_rk*mfsq)/s))))/(2.0_rk*pi*sqrt((masq-s)*(masq-s))*sqrt(-s*s*(-4.0_rk mfsq+s)*(-4.0_rk*mfsq+s)))
    !sigma_agff = gaff*gaff*nc*2.0_rk*mfsq*( mfsq* sqrt(s* (-4.0_rk*mfsq+s))*(-2.0_rk*masq*sqrt(1.0_rk-(4.0_rk* mfsq)/s)*s+(masq*masq-4.0_rk*masq*mfsq+s*s)* &
    !              log((-2.0_rk*mfsq+s+sqrt(1.0_rk-(4.0_rk*mfsq)/s)*s)/(2.0_rk*mfsq))))/(4.0_rk*pi*sqrt(1.0_rk-(4.0_rk*mfsq)/s)*((masq-s)*(masq-s))**(1.5)*s)
    sigma_agff = (gaff*gaff*mfsq*Sqrt(s*(-4.0_rk*mfsq + s))*(-2.0_rk*masq*Sqrt(1.0_rk - (4.0_rk*mfsq)/s)*s + (masq*masq - 4.0_rk*masq*mfsq + s*s)*&
                Log((-2.0_rk*mfsq + s + Sqrt(1.0_rk - (4.0_rk*mfsq)/s)*s)/(2.0_rk*mfsq))))/(8.0_rk*Pi*Sqrt(1.0_rk - (4.0_rk*mfsq)/s)*((masq - s)**2)**1.5*s)
    return
  end function sigma_agff

  real(kind=rk) function sigma_ffag(s, mfi,nc, ma, gaff)
    implicit none
    real(kind=rk), intent(in)         :: s, mfi, ma, gaff, nc
    real(kind=rk)                     :: masq, mfsq

    masq = ma*ma
    mfsq = mfi*mfi
    sigma_ffag = gaff*gaff/nc*2.0_rk*( mfsq* sqrt(s* (-4.0_rk*mfsq+s))*(-2.0_rk*masq*sqrt(1.0_rk-(4.0_rk* mfsq)/s)*s+(masq*masq-4.0_rk*masq*mfsq+s*s)* &
                  log((-2.0_rk*mfsq+s+sqrt(1.0_rk-(4.0_rk*mfsq)/s)*s)/(2.0_rk*mfsq))))/(4.0_rk*pi*sqrt(1.0_rk-(4.0_rk*mfsq)/s)*((masq-s)*(masq-s))**(1.5)*s)/(-4.0_rk*mfsq+s)*(-4.0_rk*masq+s)
    !sigma_agff = 2.0_rk*( mfsq* sqrt(s* (-4.0_rk*mfsq+s))*(-2.0_rk*masq*sqrt(1.0_rk-(4.0_rk* mfsq)/s)*s+(masq*masq-4.0_rk*masq*mfsq+s*s)* &
    !              log((-2.0_rk*mfsq+s+sqrt(1.0_rk-(4.0_rk*mfsq)/s)*s)/(2.0_rk*mfsq))))/(4.0_rk*pi*sqrt(1.0_rk-(4.0_rk*mfsq)/s)*((masq-s)*(masq-s))**(1.5)*s)
    return
  end function sigma_ffag

  real(kind=rk) function sigma_afgf(s, mf, ma, gaff)
    implicit none
    real(kind=rk), intent(in)         :: s, mf, ma, gaff
    real(kind=rk)                     :: masq, mfsq

    masq = ma*ma
    mfsq = mf*mf
    sigma_afgf = gaff*gaff*2.0_rk*(mfsq*(mfsq*mfsq-2.0_rk*mfsq*s+s*s)*sqrt(masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s))* &
                ((-mfsq+s)*(-mfsq+s)*(-mfsq+3.0_rk*s)+masq*(mfsq*mfsq-2.0_rk*mfsq*s-7.0_rk*s*s)-(2.0_rk*s*s* &
                (2.0_rk*masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s))* &
                (-log(-masq+mfsq+s-sqrt(masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s)))&
                +log(-masq+mfsq+s+sqrt(masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s)))))&
                /sqrt(masq*masq+(-mfsq+s)*(-mfsq+s)-2.0_rk*masq*(mfsq+s))))/&
                (16.0_rk*pi*(mfsq-s)*(mfsq-s)*(mfsq-s)*s*s*&
                (masq*masq-2.0_rk*masq*mfsq+mfsq*mfsq-2.0_rk*masq*s-2.0_rk*mfsq*s+s*s))
    return
  end function sigma_afgf

  real(kind=rk) function sigma_agff_th(s, mfi,nc, ma, mg, gaff)
    implicit none
    real(kind=rk), intent(in)         :: s, mfi, ma, mg, gaff, nc
    real(kind=rk)                     :: masq, mfsq, mgsq

    masq = ma*ma
    mfsq = mfi*mfi
    mgsq = mg*mg
!    sigma_agff_th = 8.0_rk*pi*gaff*gaff*((Sqrt(s*(-4.0_rk*mfsq + s))* &
!                  (-((masq*(mgsq +2.0_rk* mfsq)*(-s**(1.5_rk)*Sqrt(-4.0_rk* mfsq + s) + (masq + mgsq)*&
!                  Sqrt(s*(-4.0_rk*mfsq +s))))/((masq - mgsq)*(masq - mgsq)*mfsq + masq*mgsq*s -&
!                  2.0_rk* (masq + mgsq)*mfsq*s + mfsq*s*s)) - ((masq*masq -&
!                  4.0_rk* masq*mfsq + (mgsq - s)*(mgsq - s))*Log(((masq + mgsq - s)*&
!                  Sqrt(s) - Sqrt(-4.0_rk* mfsq + s)*Sqrt(masq*masq + (mgsq - s)*(mgsq - s) -&
!                  2.0_rk* masq*(mgsq + s)))/((masq + mgsq - s)*Sqrt(s) +&
!                  Sqrt(-4.0_rk* mfsq + s)*Sqrt(masq*masq + (mgsq - s)*(mgsq - s) - 2.0_rk* masq*(mgsq + s)))))&
!                  /Sqrt(masq*masq + (mgsq - s)*(mgsq - s) - 2.0_rk* masq*(mgsq + s))))/(16.0_rk *pi*pi*&
!                  Sqrt(1.0_rk - (4.0_rk* mfsq)/s)*(masq + mgsq - s)*s*&
!                  Sqrt(masq*masq + (mgsq - s)*(mgsq - s) - 2.0_rk* masq*(mgsq + s))))

  sigma_agff_th = -(gaff*gaff*mfsq*nc*((masq*(2.0_rk*mfsq + mgsq)*&
                  (masq + mgsq - s)*Sqrt(s*(-4.0_rk*mfsq + s)))/&
                  ((mfsq*(masq - mgsq)*(masq - mgsq) + masq*mgsq*s - &
                  2.0_rk*mfsq*(masq + mgsq)*s + mfsq*s*s)*Sqrt(masq*masq &
                  + (mgsq - s)*(mgsq - s) - 2.0_rk*masq*(mgsq + s)))&
                   + ((masq*masq - 4.0_rk*masq*mfsq + (mgsq - s)*(mgsq - s))&
                   *Log((1.0_rk + (Sqrt(1.0_rk - (4.0_rk*mfsq)/s)*&
                   Sqrt(masq*masq + (mgsq - s)*(mgsq - s)- 2*masq*(mgsq + s)))/&
                   (-masq - mgsq + s))/(1.0_rk - (Sqrt(1.0_rk - (4.0_rk*mfsq)/s)&
                   *Sqrt(masq*masq + (mgsq - s)*(mgsq - s) - &
                   2.0_rk*masq*(mgsq + s)))/(-masq - mgsq + s))))/(masq*masq + &
                   (mgsq - s)*(mgsq - s) - 2.0_rk*masq*(mgsq + s))))/(masq + mgsq - s)
                  !*alpha*q_f*qf*nc
    return
  end function sigma_agff_th

  real(kind=rk) function sigma_ahff(s, mfi, nc, ma, gaff)
    implicit none
    real(kind=rk), intent(in)   :: s, mfi, ma, gaff, nc
    real(kind=rk)               :: masq, mhsq

    mhsq = mh*mh
    masq = ma*ma
    sigma_ahff = (gaff*gaff*mfi*mfi*nc*Sqrt((-4.0_rk*mfi*mfi*s + s*s)/&
                  (masq*masq - 2.0_rk*masq*mhsq + mhsq*mhsq - 2.0_rk*masq*s&
                   - 2.0_rk*mhsq*s + s*s)))/(8.0_rk*pi*v*v)

    return
  end function sigma_ahff


!  real(kind=rk) function sigma_afgf_th(s, mfi, ma, mg, gaff, argsint)
!    implicit none
!    real(kind=rk), intent(in)         :: s, mfi, ma, mg, gaff
!    type (type_argsint), intent(inout)  :: argsint
!    real(kind=rk)                     :: masq, mfsq, mgsq
!
!    masq = ma*ma
!    mfsq = mfi*mfi
!    mgsq = mg*mg
!    sigma_afgf_th = -2.0_rk*4.0_rk*pi*mfsq*gaff*gaff*(1.0_rk/&
!                    (64.0_rk* pi*pi* (mfsq - s)*(mfsq - s))*&
!                    ((Sqrt(mgsq*mgsq + (mfsq - s)*(mfsq - s) - 2.0_rk* mgsq*(mfsq + s))&
!                    *(mfsq*(mfsq - s)**4*(mgsq - mfsq + 3.0_rk* s) -&
!                    masq*masq*masq* mgsq*(mfsq + s)*(mgsq - mfsq + 3.0_rk* s) +&
!                    masq* (mfsq - s)*(mfsq - s)*(mgsq*mgsq*mgsq + mfsq*mfsq*mfsq - &
!                    2.0_rk*mfsq*mfsq*s - 7.0_rk*mfsq*s*s +&
!                    2.0_rk* mgsq*mgsq*(-2.0_rk* mfsq + s) + &
!                    mgsq*(mfsq - 5.0_rk*s)*(2.0_rk* mfsq + s)) -&
!                    masq*masq*mgsq*(mgsq - mfsq + 3.0_rk* s)*(mgsq*(mfsq + s) -&
!                    2.0_rk* (2.0_rk* mfsq*mfsq + mfsq*s + s*s))))/&
!                    (s*s*Sqrt(masq*masq + (mfsq - s)*(mfsq - s) -&
!                    2.0_rk* masq*(mfsq + s))*(masq*masq*mgsq +&
!                    masq*mgsq*(mgsq - 3.0_rk* mfsq - s) + &
!                    (mfsq*mfi - mfi*s)*(mfsq*mfi - mfi*s))) +&
!                    1.0_rk/(masq*masq + (mfsq - s)*(mfsq - s) &
!                    - 2.0_rk* masq*(mfsq + s))*&
!                    2.0_rk* (mfsq - s)*(2.0_rk* masq*masq + (mfsq - s)*(mfsq - s) -&
!                    2.0_rk* masq*(mfsq + s))*Log(((ma - mfi)*(mg - mfi)*(ma + mfi)*(mg + mfi) + (masq +&
!                    mgsq)*s - s*s - Sqrt(masq*masq + (mfsq - s)*(mfsq - s)&
!                     - 2.0_rk* masq*(mfsq + s))*Sqrt(&
!                    mgsq*mgsq + (mfsq - s)*(mfsq - s) - 2.0_rk* mgsq*(mfsq + s)))/&
!                    ((ma - mfi)*(mg - mfi)*(ma + mfi)*(mg + mfi) + (masq + mgsq)*s - s*s +&
!                    Sqrt(masq*masq + (mfsq - s)*(mfsq - s) - 2.0_rk* masq*(mfsq + s))*Sqrt(&
!                    mgsq*mgsq + (mfsq - s)*(mfsq - s) - 2.0_rk* mgsq*(mfsq + s))))))
!
!    return
!  end function sigma_afgf_th

real(kind=rk) function sigma_afgf_th(s, mfi, ma, mg, gaff, argsint)
  implicit none
  real(kind=rk), intent(in)           :: s, mfi, ma, mg, gaff
  type (type_argsint), intent(inout)  :: argsint
  real(kind=rk)                       :: masq, mfsq, mgsq, epsabs, epsrel, result, abserr
  integer(kind=ik)                    :: neval, ier

  masq = ma*ma
  mfsq = mfi*mfi
  mgsq = mg*mg
  epsabs = 1e-20_rk
  epsrel = 1e-20_rk
  argsint%s = s
  call qags(M2_afgf_th,argsint,-1.0_rk,&
            1.0_rk, epsabs, epsrel, result, abserr, neval, ier)
  !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
  sigma_afgf_th = argsint%g*argsint%g*sqrt(lambda(s,mfsq,mgsq)/lambda(s,mfsq,masq))/64.0_rk/pi/pi/s*result
  return
end function sigma_afgf_th

real(kind=rk) function M2_afgf_th(ctheta, argsint)
  implicit none
  real(kind=rk), intent(in)           :: ctheta
  type (type_argsint), intent(in)     :: argsint
  real(kind=rk)                       :: masq, mfsq, mgsq
  real(kind=rk)                       :: mg, ma, mf, s, ssq, cthetasq, mfsqmsq
  mg = argsint%mg
  ma = argsint%ma
  mf = argsint%mf
  s = argsint%s

  mgsq = mg*mg
  masq = ma*ma
  ssq = s*s
  mfsq = mf*mf
  cthetasq = ctheta*ctheta
  mfsqmsq = (mfsq - s)*(mfsq - s)

  M2_afgf_th = -mfsq*2.0_rk*4.0_rk*pi*(masq*masq*masq*(-(mgsq*mgsq*mgsq*((1.0_rk+ 3.0_rk*cthetasq)*mfsq + s - cthetasq*s)) + &
                mgsq*mgsq*((3 + 9*cthetasq)*mfsq*mfsq - 2.0_rk*(1.0_rk+ cthetasq)*mfsq*s - &
                (5.0_rk + 3.0_rk*cthetasq)*ssq) + &
                mgsq*(-3.0_rk*(1.0_rk+ 3.0_rk*cthetasq)*mfsq*mfsq*mfsq + &
                (7.0_rk + 13.0_rk*cthetasq)*mfsq*mfsq*s + (-5.0_rk + 9.0_rk*cthetasq)*mfsq*ssq + &
                (1.0_rk+ 3.0_rk*cthetasq)*ssq*s) + mfsqmsq*((1.0_rk+ 3.0_rk*cthetasq)*mfsq*mfsq &
                - (-5 + cthetasq)*ssq - 2*mfsq*(s + 3.0_rk*cthetasq*s))) - &
                masq*(mfsq - s)*(mgsq*mgsq*mgsq*(mfsq - s)*((3.0_rk + 9.0_rk*cthetasq)*mfsq + &
                (-1.0_rk+ 5*cthetasq)*s) + mgsq*mgsq*(-9.0_rk*(1.0_rk+ 3.0_rk*cthetasq)*mfsq*mfsq*mfsq + &
                (19 + cthetasq)*mfsq*mfsq*s + (5.0_rk + 7.0_rk*cthetasq)*ssq*s - &
                2.0_rk*ctheta*(-1.0_rk+ cthetasq)*s*Sqrt(masq*masq + mfsqmsq - &
                2.0_rk*masq*(mfsq + s))*Sqrt(mgsq*mgsq + mfsqmsq - &
                2.0_rk*mgsq*(mfsq + s)) + mfsq*((-15.0_rk + 19.0_rk*cthetasq)*ssq - &
                2.0_rk*ctheta*(3.0_rk + cthetasq)*Sqrt(masq*masq + mfsqmsq - &
                2.0_rk*masq*(mfsq + s))*Sqrt(mgsq*mgsq + mfsqmsq - &
                2.0_rk*mgsq*(mfsq + s)))) + mgsq*(9.0_rk*(1.0_rk+ 3.0_rk*cthetasq)*mfsq*mfsq*mfsq*mfsq - &
                26.0_rk*(1.0_rk+ cthetasq)*mfsq*mfsq*mfsq*s + (-9.0_rk + cthetasq)*ssq*ssq + &
                4.0_rk*ctheta*(1.0_rk+ cthetasq)*ssq*Sqrt(masq*masq + mfsqmsq - &
                2.0_rk*masq*(mfsq + s))*Sqrt(mgsq*mgsq + mfsqmsq - &
                2.0_rk*mgsq*(mfsq + s)) + 2.0_rk*mfsq*s*((-3.0_rk + 5.0_rk*cthetasq)*ssq + &
                4.0_rk*ctheta*(-2.0_rk + cthetasq)*Sqrt(masq*masq + mfsqmsq - &
                2.0_rk*masq*(mfsq + s))*Sqrt(mgsq*mgsq + mfsqmsq - &
                2.0_rk*mgsq*(mfsq + s))) + 4.0_rk*mfsq*mfsq*((8.0_rk - 3.0_rk*cthetasq)*ssq + &
                ctheta*(3 + cthetasq)*Sqrt(masq*masq + mfsqmsq - &
                2.0_rk*masq*(mfsq + s))*Sqrt(mgsq*mgsq + mfsqmsq - &
                2.0_rk*mgsq*(mfsq + s)))) - mfsqmsq*((3.0_rk + 9.0_rk*cthetasq)*mfsq*mfsq*mfsq &
                - (5.0_rk + 11.0_rk*cthetasq)*mfsq*mfsq*s + (-1.0_rk+ cthetasq)*s* &
                (3.0_rk*ssq + 2.0_rk*ctheta*Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s)) &
                *Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))) + &
                 mfsq*((5.0_rk - cthetasq)*ssq + 2.0_rk*ctheta*(3.0_rk + cthetasq)* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))))) - &
                 (mfsq - s)**3*((-1.0_rk- 3.0_rk*cthetasq)*mg**6*(mfsq - s) + &
                 mgsq*mgsq*((3.0_rk + 9.0_rk*cthetasq)*mfsq*mfsq - 4*(1.0_rk+ cthetasq)*mfsq*s + &
                 (1.0_rk- 5.0_rk*cthetasq)*ssq + ctheta*(3.0_rk + cthetasq)* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))) + &
                 (mfsq - s)*((1.0_rk+ 3.0_rk*cthetasq)*mfsq*mfsq*mfsq - mfsq*mfsq*(s + 7.0_rk*cthetasq*s) &
                 - (-1.0_rk+ cthetasq)*s*(ssq + ctheta* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))) + &
                 mfsq*((-1.0_rk+ 5.0_rk*cthetasq)*ssq + ctheta*(3.0_rk + cthetasq)* &
                 Sqrt(masq*masq + mfsqmsq - 2*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2*mgsq*(mfsq + s)))) &
                 - mgsq*((3 + 9*cthetasq)*mf**6 - (5.0_rk + 11.0_rk*cthetasq)*mfsq*mfsq*s + &
                 (-1.0_rk+ cthetasq)*s*(-ssq + 2.0_rk*ctheta* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))) + &
                 mfsq*((1.0_rk+ 3.0_rk*cthetasq)*ssq + 2.0_rk*ctheta*(3.0_rk + cthetasq)* &
                 Sqrt(masq*masq + mfsqmsq - 2*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2*mgsq*(mfsq + s))))) + &
                 masq*masq*(mgsq*mgsq*mgsq*((3.0_rk + 9.0_rk*cthetasq)* &
                 mfsq*mfsq - 2.0_rk*(1.0_rk+ cthetasq)*mfsq*s + &
                 (-1.0_rk+ cthetasq)*ssq) - mgsq*mgsq* &
                 (9.0_rk*(1.0_rk+ 3.0_rk*cthetasq)*mfsq*mfsq*mfsq + &
                 (-17.0_rk + cthetasq)*mfsq*mfsq*s - (-1.0_rk+ cthetasq)*s* &
                 (ssq + ctheta*Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))) + &
                 mfsq*((7.0_rk - 11.0_rk*cthetasq)*ssq + ctheta*(3.0_rk + cthetasq)* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s)))) - &
                 mfsqmsq*((3.0_rk + 9.0_rk*cthetasq)*mf**6 - &
                 (7.0_rk + 13.0_rk*cthetasq)*mfsq*mfsq*s + &
                 mfsq*((13.0_rk - 9.0_rk*cthetasq)*ssq + &
                 ctheta*(3.0_rk + cthetasq)*Sqrt(masq*masq + &
                 mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))) - &
                 s*((-7.0_rk + 3.0_rk*cthetasq)*ssq + ctheta*(3.0_rk + cthetasq)* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s)))) + &
                 mgsq*(9.0_rk*(1.0_rk+ 3.0_rk*cthetasq)*mfsq*mfsq*mfsq*mfsq - &
                 28.0_rk*(1.0_rk+ cthetasq)*mf**6*s + &
                 mfsq*mfsq*((38.0_rk - 30.0_rk*cthetasq)*ssq &
                 + 2.0_rk*ctheta*(3.0_rk + cthetasq)* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))) - &
                 ssq*((-9.0_rk + 5.0_rk*cthetasq)*ssq + 2.0_rk*ctheta*(3.0_rk + cthetasq)* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s))) - &
                 4.0_rk*mfsq*(7.0_rk*(1.0_rk+ cthetasq)*ssq*s + 2.0_rk*ctheta*s* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s)))))) &
                 /(2.*mfsqmsq*s*((ma - mf)*(ma + mf)*(-mf + mg)* &
                 (mf + mg) + (masq + mgsq)*s - ssq + ctheta* &
                 Sqrt(masq*masq + mfsqmsq - 2.0_rk*masq*(mfsq + s))* &
                 Sqrt(mgsq*mgsq + mfsqmsq - 2.0_rk*mgsq*(mfsq + s)))**2)
  return
end function M2_afgf_th

  real(kind=rk) function lambda(x,y,z)
    implicit none
    real(kind=rk), intent(in)   :: x, y, z
    lambda = x*x + y*y + z*z - 2.0_rk*x*y - 2.0_rk*x*z - 2.0_rk*y*z
    return
  end function lambda

  real(kind=rk) function M2ffa(mfi, ma)
    implicit none
    real(kind=rk)           :: mfi, ma

    M2ffa =  ma*ma*2.0_rk
    return
  end function M2ffa

  real(kind=rk) function Gamma_ffa(mfi, ma)
    implicit none
    real(kind=rk)           :: mfi, ma

    Gamma_ffa = ma/(8.0_rk*pi)*sqrt(1.0_rk-4.0_rk*mfi*mfi/ma/ma)
    return
  end function Gamma_ffa


!  real(kind=rk) function Mgga()
!    implicit none
!
!    Mgga = 2.0_rk/pi/pi*ma*ma*ma*ma*
!    Total[mFin[[1, 1, All]]^2 ncf[[1, All]] qf^2 Table[
!2 Abs[ArcSin[Sqrt[ma^2/(4*mFin[[1, 1, i]]^2)]]^2]/ma^2, {i, 1,
!Length[mFin[[1, 1, All]]]}]]^2
!
!  end function Mgga
end module
