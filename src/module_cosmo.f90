module module_cosmo
  use module_precision
  use module_utils
  use module_params
  use module_dof

  contains

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
      call interp_linear(nd, params%geff_HS(1,:),params%geff_HS(2,:),T, gHS)
      call geffSM(T,params,gSM)
      ! Axion temperature
      nr = size(params%rhoa_rho,2)
      call interp_linear(nr, params%rhoa_rho(1,:),params%rhoa_rho(2,:),T, rar)
      !Ta = sqrt(sqrt(( gSM*gSM/gHS *rar )))*T
      Ta = sqrt(sqrt(( gSM*gSM/(ga) *rar )))*T
    end function Ta


    include "initial_conditions.f90"
end module
