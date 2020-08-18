module module_xsecs
use module_precision
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
end module
