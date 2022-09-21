! effective degrees of freedom (aso for the HS). SM effective dofs are in geff_new.f90
module module_dof
  use module_precision
  use module_utils
  use module_params

  contains
    real(kind=rk) function heff(T,params)

      implicit none
      type (type_params), intent(in)        :: params
      real(kind=rk), intent(in)             :: T
      real(kind=rk)                         :: hSM, hHS
      integer(kind=ik)                      :: nd

      nd = size(params%heff_HS,2)
      !call heffSM(T,params,gSM)
      hSM = geff_s(T)
      call interp_linear(nd, params%heff_HS(1,:),params%heff_HS(2,:),log10(T), hHS)
      heff = hSM + hHS
      return
    end function heff

    real(kind=rk) function geff(T,params)

      implicit none
      type (type_params), intent(in)        :: params
      real(kind=rk), intent(in)             :: T
      real(kind=rk)                         :: gSM, gHS
      integer(kind=ik)                      :: nd

      nd = size(params%geff_HS,2)
      gSM = geff_rho(T)
      call interp_linear(nd, params%geff_HS(1,:),params%geff_HS(2,:),log10(T), gHS)
      geff = gSM + gHS
      return
    end function geff

    real(kind=rk) function sqrtgstar(T,params)

      implicit none
      type (type_params), intent(in)        :: params
      real(kind=rk), intent(in)             :: T
      real(kind=rk)                         :: gSM, gHS, hHS, hSM
      integer(kind=ik)                      :: nd,nd2

      nd = size(params%geff_HS,2)
      nd2 = size(params%heff_HS,2)
      gSM = geff_rho(T)
      hSM = geff_s(T)
      call interp_linear(nd, params%geff_HS(1,:),params%geff_HS(2,:),log10(T), gHS)
      call interp_linear(nd2, params%heff_HS(1,:),params%heff_HS(2,:),log10(T), hHS)
      sqrtgstar = sqrt(gSM) *(1 + hHS / hSM - 0.5_rk*gSM*gHS / hSM/hSM)
      return
    end function sqrtgstar

    real(kind=rk) function hi(y,argsint)
      implicit none
        real(kind=rk), intent(in)           :: y
        type (type_argsint), intent(in)     :: argsint
        real(kind=rk)                       :: x
        x = argsint%mf/argsint%T
        hi = 45.0_rk*argsint%g/4.0_rk/pi/pi/pi/pi* x*x*x*x*&
                    y * sqrt( y*y - 1.0_rk )/(exp(x*y)+dble((-1)**mod(int(argsint%g),2))) * (4.0_rk*y*y-1.0_rk)/3.0_rk/y

    end function hi

    real(kind=rk) function gi(y,argsint)
      implicit none
        real(kind=rk), intent(in)           :: y
        type (type_argsint), intent(in)     :: argsint
        real(kind=rk)                       :: x
        x = argsint%mf/argsint%T
        gi = 15.0_rk*argsint%g/pi/pi/pi/pi* x*x*x*x*&
            y * sqrt( y*y - 1.0_rk )/(exp(x*y)+dble((-1)**mod(int(argsint%g),2))) * y

    end function gi

    subroutine geffHS(params,argsint,geff)
      implicit none
      type (type_params), intent(in)              :: params
      type (type_argsint), intent(inout)          :: argsint
      real(kind=rk), dimension(:,:), intent(out)  :: geff
      integer(kind=ik)                            :: k, ier, neval
      real(kind=rk)                               :: result, T, abserr, epsabs, epsrel, facTdi

      ! factor for Tdi: Tdi = m/facTdi
      facTdi = 20.0_rk
      geff = 0.0_rk
      epsabs=1e-8_rk
      epsrel=1e-5_rk
      call linspace(-2.0_rk,2.0_rk,geff(1,:))
      do k=1,size(geff,2)
        T = 10.0_rk**geff(1,k)
        !geff(1,k) = T
        argsint%mf = params%ma
        argsint%g = ga
        if (T>=params%ma/facTdi) then
          argsint%T = T
          call qagi(gi,argsint,1.0_rk,&
                  1, epsabs, epsrel, result, abserr, neval, ier)
        else
          argsint%T = params%ma/facTdi
          call qagi(gi,argsint,1.0_rk,&
                  1, epsabs, epsrel, result, abserr, neval, ier)
        end if
        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr,T
        geff(2,k) = geff(2,k) + result
        argsint%mf = params%mx
        argsint%g = gDM
        if (T>=params%mx/facTdi) then
          argsint%T = T
          call qagi(gi,argsint,1.0_rk,&
                  1, epsabs, epsrel, result, abserr, neval, ier)
        else
          argsint%T = params%mx/facTdi
          call qagi(gi,argsint,1.0_rk,&
                  1, epsabs, epsrel, result, abserr, neval, ier)
        end if
        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr,T
        geff(2,k) = geff(2,k) + result
      end do

    end subroutine geffHS


!    subroutine geffHS(params,argsint,geff)
!      implicit none
!      type (type_params), intent(in)              :: params
!      type (type_argsint), intent(inout)          :: argsint
!      real(kind=rk), dimension(:,:), intent(out)  :: geff
!      integer(kind=ik)                            :: k, ier, neval
!      real(kind=rk)                               :: result, T, abserr, epsabs, epsrel
!      real(kind=rk)                               :: heff, heffTdx, hx, hxTdx, gx, facTdi
!
!      facTdi = 20.0_rk
!      geff = 0.0_rk
!      epsabs=1e-8_rk
!      epsrel=1e-5_rk
!      call linspace(-2.0_rk,2.0_rk,geff(1,:))
!      do k=1,size(geff,2)
!        T = 10.0_rk**geff(1,k)
!        !geff(1,k) = T
!        argsint%mf = params%ma
!        argsint%g = ga
!        if (T>=params%ma/facTdi) then
!          argsint%T = T
!          call qagi(gi,argsint,1.0_rk,&
!                  1, epsabs, epsrel, result, abserr, neval, ier)
!        else
!          argsint%T = params%ma/facTdi
!          call qagi(gi,argsint,1.0_rk,&
!                  1, epsabs, epsrel, result, abserr, neval, ier)
!        end if
!        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr,T
!        geff(2,k) = result
!
!        argsint%mf = params%mx
!        argsint%g = gDM
!        if (T>=params%mx/facTdi) then
!          argsint%T = T
!          call qagi(gi,argsint,1.0_rk,&
!                  1, epsabs, epsrel, result, abserr, neval, ier)
!          if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr,T
!          geff(2,k) = geff(2,k) + result
!        else
!          call qagi(gi,argsint,1.0_rk,&
!                  1, epsabs, epsrel, result, abserr, neval, ier)
!          if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr,T
!          gx = result
!          call qagi(hi,argsint,1.0_rk,&
!                  1, epsabs, epsrel, result, abserr, neval, ier)
!          if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr,T
!          hx = result
!          call interp_linear(size(params%heff_HS,2), params%heff_HS(1,:),params%heff_HS(2,:),geff(1,k),heff)!T, heff)
!          call interp_linear(size(params%heff_HS,2), params%heff_HS(1,:),params%heff_HS(2,:),log10(params%mx/facTdi), heffTdx)
!          argsint%T = params%mx/facTdi
!          call qagi(hi,argsint,1.0_rk,&
!                  1, epsabs, epsrel, result, abserr, neval, ier)
!          if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr,T
!          hxTdx = result
!          geff(2,k) = geff(2,k) + gx * (heff*hxTdx/hx/heffTdx)**(4.0_rk/3.0_rk)
!        end if
!      end do
!
!    end subroutine geffHS


    subroutine heffHS(params,argsint,heff)
      implicit none
      type (type_params), intent(in)              :: params
      type (type_argsint), intent(inout)          :: argsint
      real(kind=rk), dimension(:,:), intent(out)  :: heff
      integer(kind=ik)                            :: i, ier, neval, k
      real(kind=rk)                               :: result, T, abserr, facTdi
      real(kind=rk)                               :: epsabs, epsrel, hc, hxTi, haTi, hcTa, hcTx

      facTdi = 20.0_rk
      heff = 0.0_rk
      epsabs=1e-8_rk
      epsrel=1e-5_rk
      call linspace(-2.0_rk,2.0_rk,heff(1,:))
      do k=1,size(heff(1,:))
        hc = 0.0_rk
        hxTi = 0.0_rk
        haTi = 0.0_rk
        hcTx = 0.0_rk
        hcTa = 0.0_rk
        T = 10.0_rk**heff(1,k)
        !heff(1,k) = T
        if (T<params%mx/facTdi .and. T<params%ma/facTdi) then
          heff(2,k) = 0.0_rk
        else
          argsint%mf = params%ma
          argsint%g = ga
          if (T>=params%ma/facTdi) then
            argsint%T = T
            call qagi(hi,argsint,1.0_rk,&
                    1, epsabs, epsrel, result, abserr, neval, ier)
            hc = hc + result
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr
            hcTa = 1.0_rk
          else
            argsint%T = params%ma/facTdi
            call qagi(hi,argsint,1.0_rk,&
                    1, epsabs, epsrel, result, abserr, neval, ier)
            haTi = result
            hcTa = hcTa + result
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr
            argsint%mf = params%mx
            argsint%g = gDM
            call qagi(hi,argsint,1.0_rk,&
                    1, epsabs, epsrel, result, abserr, neval, ier)
            hcTa = hcTa + result
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr
          end if
          argsint%mf = params%mx
          argsint%g = gDM
          if (T>=params%mx/facTdi) then
            argsint%T = T
            call qagi(hi,argsint,1.0_rk,&
                    1, epsabs, epsrel, result, abserr, neval, ier)
            hc = hc + result
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr
            hcTx = 1.0_rk
          else
            argsint%T = params%mx/facTdi
            call qagi(hi,argsint,1.0_rk,&
                    1, epsabs, epsrel, result, abserr, neval, ier)
            hxTi = result
            hcTx = hcTx + result
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr
            argsint%mf = params%ma
            argsint%g = ga
            call qagi(hi,argsint,1.0_rk,&
                    1, epsabs, epsrel, result, abserr, neval, ier)
            hcTx = hcTx + result
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr
          end if
          heff(2,k) = hc * (1.0_rk + haTi/hcTa) * (1.0_rk +  hxTi/hcTx)
        end if
      end do

    end subroutine heffHS

    include "ini_cons_to_params.f90"
    include "geff_new.f90"

  end module module_dof
