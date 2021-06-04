module module_rhs

  use mpi
  use module_precision
  use module_params
  use module_utils
  use module_cosmo
  use module_xsecs

  contains

    real(kind=rk) function gammav(T, params, Gamma)
      implicit none
      type (type_params), intent(in)        :: params
      real(kind=rk), intent(in)             :: T
      character(len=*), intent(in)          :: Gamma
      real(kind=rk)                         :: mfi
      integer(kind=ik)                      :: i

      gammav = 0.0_rk
      select case(Gamma)
        case("affth")
          do i=1, 3
            mfi = ml_th(T,mf(i),qf(i))
            if (2.0_rk*mfi<params%ma) then
              gammav = gammav + Gamma_ffa(mfi,params%ma)
            end if
          end do
          if (T>QCDcut) then
            do i=4,9
              mfi = mq_th(T,mf(i),qf(i))
              if (2.0_rk*mfi<params%ma) then
                gammav = gammav + Gamma_ffa(mfi,params%ma)
              end if
            end do
          end if
        case("aff")
          do i=1, 3
            if (2.0_rk*mf(i)<params%ma) then
              gammav = gammav + Gamma_ffa(mf(i),params%ma)
            end if
          end do
          if (T>QCDcut) then
            do i=4,9
              if (2.0_rk*mf(i)<params%ma) then
                gammav = gammav + Gamma_ffa(mf(i),params%ma)
              end if
            end do
          end if
        case default
          write(*,*) "Error! decay width ", Gamma, " not implemented"
        end select

        if (T/params%ma>0.03_rk) then
          gammav = gammav*bessK1(params%ma/T)/bessK2(params%ma/T)
        else
          gammav = gammav*(1.0_rk-3.0_rk/2.0_rk*T/params%ma+15.0_rk*T*T/8.0_rk/params%ma/params%ma)
        end if
    end function gammav
    subroutine gamma_r_new( T, params, argsint, sigma, gam )
      implicit none
      type (type_params), intent(in)        :: params
      type (type_argsint), intent(inout)    :: argsint
      real(kind=rk), intent(in)             :: T
      character(len=*), intent(in)          :: sigma
      real(kind=rk), intent(out)            :: gam
      real(kind=rk)                         :: mx, result, epsabs, epsrel, abserr,ma, alphas, alphaem
      integer(kind=ik)                      :: i, ier, neval, nd

      epsabs = 1e-7_rk
      epsrel = 1e-2_rk
      mx = params%mx
      ma = params%ma
      gam = 0.0_rk
      if (T>QCDcut) then
        alphas = alpha_s(2.0_rk*pi*T)
      end if
      alphaem = alpha_qed_th(2.0_rk*pi*T)
      !end if
      argsint%T=T

      select case(sigma)
      case("agff")
        do i=1, 3
          argsint%mf = mf(i)
          call qagi(kernel_agff,argsint,max(ma*ma,4.0_rk*mf(i)*mf(i)),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            !write(*,*) "Integral did not converge for the first tim ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            call qags(kernel_agff,argsint,max(ma*ma,4.0_rk*mf(i)*mf(i)),&
                      1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
            !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + alphaem*qf(i)*qf(i)*ga*ggamma*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
        if (T>QCDcut) then
          do i=4, 9
            argsint%mf = mf(i)
            call qagi(kernel_agff,argsint,max(QCDcut*QCDcut,max(ma*ma,4.0_rk*mf(i)*mf(i))),&
                      1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_agff,argsint,max(QCDcut*QCDcut,max(ma*ma,4.0_rk*mf(i)*mf(i))),&
                        1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
            !if (ier > 0)  write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            end if
            gam = gam + (0.5_rk*alphas*ga*gg+alphaem*qf(i)*qf(i)*ga*ggamma)&
                        *T/(32.0_rk*pi*pi*pi*pi)*result
          end do
        end if
      case("afgf")
        do i=1, 3
          argsint%mf = mf(i)
          call qagi(kernel_afgf,argsint,(ma+mf(i))*(ma+mf(i)),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            call qags(kernel_afgf,argsint,(ma+mf(i))*(ma+mf(i)),&
                      1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + alphaem*qf(i)*qf(i)*ga*gf(i)*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
        if (T>QCDcut) then
          do i=4, 9
            argsint%mf = mf(i)
            call qagi(kernel_afgf,argsint,max(QCDcut*QCDcut,(ma+mf(i))*(ma+mf(i))),&
                      1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_afgf,argsint,max(QCDcut*QCDcut,(ma+mf(i))*(ma+mf(i))),&
                        1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
              if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            end if
            gam = gam + (4.0_rk/3.0_rk*alphas*ga*gf(i)+alphaem*qf(i)*qf(i)*ga*gf(i))&
                        *T/(32.0_rk*pi*pi*pi*pi)*result
          end do
        end if
      case("xxff")
        do i=1, 9
          if (i>3 .and. T<QCDcut) then
          else
            argsint%mf = mf(i)
            argsint%nc = ncf(i)
            call qagi(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*mf(i)*mf(i)),&
                      1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*mf(i)*mf(i)),&
                        1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
              !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            end if
            gam = gam + gDM*gDM*T/(32.0_rk*pi*pi*pi*pi)*result
          end if
        end do
      case("xxffth")
        do i=1, 3
          argsint%mf = ml_th(T,mf(i),qf(i))
          argsint%nc = ncf(i)
          call qagi(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*argsint%mf*argsint%mf),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            call qags(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*argsint%mf*argsint%mf),&
                      1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
            !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + gDM*gDM*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
        if (T>QCDcut) then
          do i=4,9
            argsint%mf = mq_th(T,mf(i),qf(i))
            argsint%nc = ncf(i)
            call qagi(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*argsint%mf*argsint%mf),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*argsint%mf*argsint%mf),&
                      1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
            !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            end if
            gam = gam + gDM*gDM*T/(32.0_rk*pi*pi*pi*pi)*result
          end do
        end if
      case("agffth")
        argsint%mg = mgamma_th(T)
        do i=1, 3
          argsint%mf = ml_th(T,mf(i),qf(i))
          call qagi(kernel_agff_th,argsint,max((ma+argsint%mg)*&
                    (ma+argsint%mg),4.0_rk*argsint%mf*argsint%mf),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            call qags(kernel_agff_th,argsint,max((ma+argsint%mg)*&
                      (ma+argsint%mg),4.0_rk*argsint%mf*argsint%mf),&
                      1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
            !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + alphaem*qf(i)*qf(i)*ga*ggamma*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
        if (T>QCDcut) then
          do i=4,9
            argsint%mf = mq_th(T,mf(i),qf(i))
            call qagi(kernel_agff_th,argsint,max((ma+argsint%mg)*&
                      (ma+argsint%mg),4.0_rk*argsint%mf*argsint%mf),&
                      1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_agff_th,argsint,max((ma+argsint%mg)*&
                        (ma+argsint%mg),4.0_rk*argsint%mf*argsint%mf),&
                        1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
              !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            end if
            gam = gam + alphaem*qf(i)*qf(i)*ga*ggamma*T/(32.0_rk*pi*pi*pi*pi)*result
          end do
          argsint%mg = mg_th(T)
          do i=4, 9
            argsint%mf = mq_th(T,mf(i),qf(i))
            call qagi(kernel_agff_th,argsint,max(QCDcut*QCDcut,max((ma+argsint%mg)*&
                      (ma+argsint%mg),4.0_rk*argsint%mf*argsint%mf)),&
                      1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_agff_th,argsint,max(QCDcut*QCDcut,max((ma+argsint%mg)*&
                        (ma+argsint%mg),4.0_rk*argsint%mf*argsint%mf)),&
                        1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
              !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            end if
            gam = gam + 0.5_rk*alphas*ga*gg*T/(32.0_rk*pi*pi*pi*pi)*result
          end do
        end if
      case("afgfth")
        argsint%mg = mgamma_th(T)
        do i=1, 3
          argsint%mf = ml_th(T,mf(i),qf(i))
          call qagi(kernel_afgf_th,argsint,max((ma+argsint%mf)*(ma+argsint%mf),&
                    (argsint%mg+argsint%mf)*(argsint%mg+argsint%mf)),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            call qags(kernel_afgf_th,argsint,max((ma+argsint%mf)*(ma+argsint%mf),&
                      (argsint%mg+argsint%mf)*(argsint%mg+argsint%mf)),&
                      1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + alphaem*qf(i)*qf(i)*ga*gf(i)*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
        if (T>QCDcut) then
          do i=4,9
            argsint%mf = mq_th(T,mf(i),qf(i))
            call qagi(kernel_afgf_th,argsint,max((ma+argsint%mf)*(ma+argsint%mf),&
                      (argsint%mg+argsint%mf)*(argsint%mg+argsint%mf)),&
                      1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_afgf_th,argsint,max((ma+argsint%mf)*(ma+argsint%mf),&
                        (argsint%mg+argsint%mf)*(argsint%mg+argsint%mf)),&
                        1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
              if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            end if
            gam = gam + alphaem*qf(i)*qf(i)*ga*gf(i)*T/(32.0_rk*pi*pi*pi*pi)*result
          end do
          argsint%mg = mg_th(T)
          do i=4, 9
            argsint%mf = mq_th(T,mf(i),qf(i))
            call qagi(kernel_afgf_th,argsint,max(QCDcut*QCDcut,max((ma+argsint%mf)*(ma+argsint%mf),&
                      (argsint%mg+argsint%mf)*(argsint%mg+argsint%mf))),&
                      1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_afgf_th,argsint,max(QCDcut*QCDcut,max((ma+argsint%mf)*(ma+argsint%mf),&
                        (argsint%mg+argsint%mf)*(argsint%mg+argsint%mf))),&
                        1e8_rk, epsabs, epsrel, result, abserr, neval, ier)
              if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            end if
            gam = gam + 4.0_rk/3.0_rk*alphas*ga*gf(i)*T/(32.0_rk*pi*pi*pi*pi)*result
          end do
        end if
!      case("aaxx")
!        call qagi(kernel_aaxx,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
!                  1, epsabs, epsrel, result, abserr, neval, ier)
!        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
!        gam = gam + gDM*gDM*T/(32.0_rk*pi*pi*pi*pi)*result
!      case("xxaa")
!        call qagi(kernel_xxaa,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
!                  1, epsabs, epsrel, result, abserr, neval, ier)
!        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
!        gam = gam + gDM*gDM*T/(32.0_rk*pi*pi*pi*pi)*result
      case default
        write(*,*) "Error! x section", sigma, "not implemented"
      end select

    end subroutine gamma_r_new

    subroutine sigmav( T, params, argsint, sigma, sv )
      implicit none
      type (type_params), intent(in)        :: params
      type (type_argsint), intent(inout)    :: argsint
      real(kind=rk), intent(in)             :: T
      character(len=*), intent(in)          :: sigma
      real(kind=rk), intent(out)            :: sv
      real(kind=rk)                         :: mx, result, epsabs, epsrel, abserr,ma
      integer(kind=ik)                      :: i, ier, neval, nd

      epsabs = 1e-30_rk
      epsrel = 1e-10_rk
      mx = params%mx
      ma = params%ma
      argsint%T=T

      select case(sigma)
        case("aaxx")
          if (T>0.03_rk) then
            call qagi(kernel_aaxx,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                  1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_aaxx,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                    1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
              !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
            end if
            sv = result/4.0_rk/(2.0_rk*T*ma*ma*ma*ma*bessK2(ma/T)* bessK2(ma/T))
          else
            call qagi(kernel_aaxx_series,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                  1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_aaxx_series,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                    1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
              !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
            end if
            sv = result
          end if
        case("xxaa")
          if (T>0.03_rk) then
            call qagi(kernel_xxaa,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                  1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_xxaa,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                    1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
              !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
            end if
            sv = result/4.0_rk/(2.0_rk*T*mx*mx*mx*mx* bessK2(mx/T)* bessK2(mx/T))
          else
            call qagi(kernel_xxaa_series,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                  1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_xxaa_series,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                    1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
              !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
            end if
            sv = result
          end if
        case default
          write(*,*) "Error! x section", sigma, "not implemented"
      end select
    end subroutine sigmav

    real(kind=rk) function kernel_xxff( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint

        kernel_xxff = sigma_xxff(s,argsint%mf,argsint%nc,argsint%mx,argsint%ma,argsint%g)*4.0_rk*&
                      F(s,argsint%mx,argsint%mx)*F(s,argsint%mx,argsint%mx)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel_xxff

    real(kind=rk) function drho_decay(T, ma, process)
      implicit none
      real(kind=rk), intent(in)               :: T, ma
      character(len=*), intent(in)            :: process
      real(kind=rk)                           :: mfi
      integer(kind=ik)                        :: i

      drho_decay = 0.0_rk
      select case(process)
      case("ffath")
        do i=1, 3
          mfi = ml_th(T,mf(i),qf(i))
          if (2.0_rk*mfi<ma) then
            drho_decay = drho_decay + gf(i)*gf(i)/(32.0_rk*pi*pi*pi) * ma * &
                        sqrt(ma*ma-4.0_rk*mfi*mfi)*M2ffa(mfi,ma)
          end if
        end do
        if (T>QCDcut) then
          do i=4, 9
            mfi = mq_th(T,mf(i),qf(i))
            if (2.0_rk*mfi<ma) then
              drho_decay = drho_decay + gf(i)*gf(i)/(32.0_rk*pi*pi*pi) * ma * &
                        sqrt(ma*ma-4.0_rk*mfi*mfi)*M2ffa(mfi,ma)
            end if
          end do
        end if
      case("ffa")
        do i=1, 9
          if (i>3 .and. T<QCDcut) then
          else
            if (2.0_rk*mf(i)<ma) then
              drho_decay = drho_decay + gf(i)*gf(i)/(32.0_rk*pi*pi*pi) * ma * &
                        sqrt(ma*ma-4.0_rk*mf(i)*mf(i))*M2ffa(mf(i),ma)
            end if
          end if
        end do
!        case("gga")
!          M2 = M2gga()
        case default
          write(*,*) "Error! process", process, "not implemented"
      end select
      drho_decay = drho_decay * bessK2(ma/T) * T
      return
    end function drho_decay

    include "rhs_contributions.f90"
    include "region3aeq.f90"
    include "region3a_in_n.f90"
    include "RK4.f90"
    include "rhs_contributions_in_n.f90"
    include "sm_alps.f90"
    include "thermal_masses.f90"
    include "HS_interaction.f90"
    include "initial_conditions.f90"
    include "rhop_over_rho.f90"
    include "region_freeze_out.f90"
    include "region_freeze_out_coupled.f90"
    include "region_freeze_in.f90"
    include "choose_regime.f90"
    include "rhs_contributions_general.f90"
    include "general_rhs.f90"

end module
