module module_cosmo
  use module_precision
  use module_utils
  use module_params

  contains
    subroutine geffSM(T,params,g)

      implicit none
      type (type_params), intent(in)      :: params
      real(kind=rk), intent(in)           :: T
      real(kind=rk), intent(out)          :: g
      integer(kind=ik)                    :: nd

      nd = size(params%geff,2)

      if (T > params%geff(1,nd)) then
        g = params%geff(1,nd)
      else
        call interp_linear(nd, params%geff(1,:),params%geff(2,:),T, g)
      end if
    end subroutine geffSM

    subroutine heffSM(T,params,g)

      implicit none
      type (type_params), intent(in)        :: params
      real(kind=rk), intent(in)             :: T
      real(kind=rk), intent(out)            :: g
      integer(kind=ik)                      :: nd

      nd = size(params%heff,2)

      if (T > params%heff(1,nd)) then
        g = params%heff(1,nd)
      else
        call interp_linear( nd, params%heff(1,:), params%heff(2,:), T, g )
      end if
    end subroutine heffSM
    real(kind=rk) function heff(T,params)

      implicit none
      type (type_params), intent(in)        :: params
      real(kind=rk), intent(in)             :: T
      real(kind=rk)                         :: gSM, gHS
      integer(kind=ik)                      :: nd

      nd = size(params%heff_HS,2)
      call heffSM(T,params,gSM)
      call interp_linear(nd, params%heff_HS(1,:),params%heff_HS(2,:),T, gHS)
      heff = gSM + gHS
      return
    end function heff

    real(kind=rk) function geff(T,params)

      implicit none
      type (type_params), intent(in)        :: params
      real(kind=rk), intent(in)             :: T
      real(kind=rk)                         :: gSM, gHS, hHS, hSM
      integer(kind=ik)                      :: nd

      nd = size(params%geff_HS,2)
      call geffSM(T,params,gSM)
      call heffSM (T,params,hSM)
      call interp_linear(nd, params%geff_HS(1,:),params%geff_HS(2,:),T, gHS)
      call interp_linear(nd, params%heff_HS(1,:),params%heff_HS(2,:),T, hHS)
      geff = gSM *(1 + hHS / hSM - 0.5_rk*gSM*gSM*gHS / hSM/hSM)

      !geff2rootSM[T](1+heffHS/gsSM-1/2 geff2rootSM^2 geffHS/gsSM[T]^2)
      return
    end function geff

    real(kind=rk) function Hub( T, params )
      implicit none
        real(kind=rk), intent(in)           :: T
        type (type_params), intent(in)      :: params

        Hub = 1.67 * geff(T,params) * T* T / params%Mpl
        return
    end function Hub

    real(kind=rk) function ent( T,params )
      implicit none
        real(kind=rk), intent(in)           :: T
        type (type_params), intent(in)      :: params

        ent = 2. * pi* pi /45. * heff(T,params) * T*T*T

    end function ent

    real(kind=rk) function neq( T, m, g )
      implicit none
      real(kind=rk), intent(in)           :: T, m, g

      neq = 0.5_rk * g * T * m*m * bessK2(m/T) /pi/pi

    end function neq
    real(kind=rk) function Yeq( T, params )
      implicit none
        real(kind=rk), intent(in)           :: T
        type (type_params), intent(in)      :: params
        real(kind=rk)                       :: z

          z = params%mx/T
          Yeq = 45. * 4. / (4. * pi*pi*pi*pi * heff(T,params) ) * z*z*bessK2(z)

    end function Yeq

    real(kind=rk) function Ta(T,params)
      implicit none
      real(kind=rk), intent(in)           :: T
      type (type_params), intent(in)      :: params
      integer(kind=ik)                    :: nd, nr
      real(kind=rk)                       :: gHS, gSM, rar

      nd = size(params%heff_HS,2)
      ! DOFS
      call interp_linear(nd, params%heff_HS(1,:),params%heff_HS(2,:),T, gHS)
      call geffSM(T,params,gSM)
      ! Axion temperature
      nr = size(params%rhoa_rho,2)
      call interp_linear(nr, params%rhoa_rho(1,:),params%rhoa_rho(2,:),T, rar)
      Ta = sqrt(sqrt(( gSM*gSM/gHS *rar )))*T

    end function Ta

    include "initial_conditions.f90"
end module
