module module_rhs

  use module_precision
  use module_params
  use module_utils
  use module_cosmo
  use module_xsecs

  contains

    real(kind=rk) function kernel( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint

        kernel = sigma_xxff(s,argsint%mf,argsint%nc,argsint%mx,argsint%ma)*&
                ( s - 4.0_rk* argsint%mx* argsint%mx ) *s *argsint%T * bessK2( sqrt(s)/argsint%T )
        return
    end function kernel

    real(kind=rk) function kernel2_xxaa( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint

        kernel2_xxaa = sigma_xxaa(s,argsint%mf,argsint%nc,argsint%mx,argsint%ma)* &
                      (s - 4.0_rk* argsint%mx * argsint%mx) * sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel2_xxaa

    real(kind=rk) function kernel2_xxff( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint

        kernel2_xxff = sigma_xxff(s,argsint%mf,argsint%nc,argsint%mx,argsint%ma)*&
                      (s - 4.0_rk* argsint%mx * argsint%mx) * sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel2_xxff

    subroutine gamma_r( T, params, argsint, con, gam )
      implicit none
      type (type_params), intent(in)        :: params
      type (type_argsint), intent(inout)    :: argsint
      real(kind=rk), intent(in)             :: T
      logical, intent(in)                   :: con
      real(kind=rk), intent(out)            :: gam
      real(kind=rk)                         :: mx, result, epsabs, epsrel, abserr
      integer(kind=ik)                      :: i, ier, neval

      epsabs = 1e-5_rk
      epsrel = 1e-5_rk
      gam    = 0.0_rk
      argsint%T = T
      mx = params%mx

      if ( con ) then
        do i=1, size(params%mf)
          argsint%mf = params%mf(i)
          argsint%nc = params%ncf(i)
          call qagi( kernel2_xxff, argsint, 4.0_rk*max(mx*mx,params%mf(i)*params%mf(i)), &
                    1, epsabs, epsrel, result, abserr, neval, ier )
          gam = gam + T/2.0_rk/pi**4 * result
        end do
      else
        argsint%mf = params%ma
        argsint%nc = 1.0_rk
        !!!!!!!!!!!!!!!! check this!!!!!!!
        call qagi( kernel2_xxaa, argsint, 4.0_rk*max(mx*mx,params%ma*params%ma), &
                  1, epsabs, epsrel, result, abserr, neval, ier )
        gam = T/2.0_rk/pi**4 * result
      end if
      return
    end subroutine gamma_r

    include "rhs_boltzmann.f90"
    include "RK4.f90"
end module
