module module_cosmo
  use module_precision
  use module_utils
  use module_params
  use module_dof

  contains

    real(kind=rk) function Hub_old( T, params )
      implicit none
        real(kind=rk), intent(in)           :: T
        type (type_params), intent(in)      :: params

        Hub_old = 1.67_rk * sqrtgstar(T,params) * T* T / Mpl
        return
    end function Hub_old

    real(kind=rk) function Hub( T, rhoHS)
      implicit none
        real(kind=rk), intent(in)           :: T, rhoHS
        real(kind=rk)                       :: rhoSM

        rhoSM = geff_rho(T)*pi*pi/30.0_rk*T*T*T*T
        Hub = sqrt(pi*(rhoSM+rhoHS)*8.0_rk/3.0_rk) / Mpl
        !Hub = 2.0_rk/3.0_rk*sqrt(pi*geff(T,params)/5.0_rk) * pi * T*T / Mpl
    end function Hub


    real(kind=rk) function ent( T, params )
      implicit none
        real(kind=rk), intent(in)           :: T
        type (type_params), intent(in)      :: params

        ent = 2.0_rk * pi* pi /45.0_rk * heff(T,params) * T*T*T
    end function ent

    real(kind=rk) function rho_SM(T)
      implicit none
      real(kind=rk), intent(in)   :: T

      rho_SM = pi*pi/30.0_rk*geff_rho(T)*T*T*T*T
    end function rho_SM

    include "maxwell_boltzmann.f90"
    include "HS_temperature.f90"
end module
