module module_rhs

  use mpi
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

        kernel = sigma_xxff(s,argsint%mf,argsint%nc,argsint%mx,argsint%ma,argsint%g)*&
                ( s - 4.0_rk* argsint%mx* argsint%mx ) *s *argsint%T * bessK2( sqrt(s)/argsint%T )
        return
    end function kernel

    real(kind=rk) function kernel2_xxaa( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint

        kernel2_xxaa = sigma_xxaa(s,argsint%mx,argsint%ma,argsint%g)* &
                      (s - 4.0_rk* argsint%mx * argsint%mx) * sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel2_xxaa

    real(kind=rk) function kernel2_xxff( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint

        kernel2_xxff = sigma_xxff(s,argsint%mf,argsint%nc,argsint%mx,argsint%ma,argsint%g)*&
                      (s - 4.0_rk* argsint%mx * argsint%mx) * sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel2_xxff

!    subroutine gamma_r( T, params, argsint, con, gam )
!      implicit none
!      type (type_params), intent(in)        :: params
!      type (type_argsint), intent(inout)    :: argsint
!      real(kind=rk), intent(in)             :: T
!      logical, intent(in)                   :: con
!      real(kind=rk), intent(out)            :: gam
!      real(kind=rk)                         :: mx, result, epsabs, epsrel, abserr
!      integer(kind=ik)                      :: i, ier, neval
!
!      epsabs = 1e-5_rk
!      epsrel = 1e-5_rk
!      gam    = 0.0_rk
!      argsint%T = T
!      mx = params%mx
!
!      if ( con ) then
!        do i=1, size(mf)
!          argsint%mf = mf(i)
!          argsint%nc = ncf(i)
!          call qagi( kernel2_xxff, argsint, 4.0_rk*max(mx*mx,mf(i)*mf(i)), &
!                    1, epsabs, epsrel, result, abserr, neval, ier )
!          gam = gam + T/2.0_rk/pi**4 * result
!        end do
!      else
!        !argsint%mf = params%ma
!        argsint%nc = 1.0_rk
!        !!!!!!!!!!!!!!!! check this!!!!!!!
!        call qagi( kernel2_xxaa, argsint, 4.0_rk*max(mx*mx,params%ma*params%ma), &
!                  1, epsabs, epsrel, result, abserr, neval, ier )
!        gam = T/2.0_rk/pi**4 * result
!      end if
!      return
!    end subroutine gamma_r

    subroutine gamma_r_new( T, params, argsint, sigma, gam )
      implicit none
      type (type_params), intent(in)        :: params
      type (type_argsint), intent(inout)    :: argsint
      real(kind=rk), intent(in)             :: T
      character(len=*), intent(in)          :: sigma
      real(kind=rk), intent(out)            :: gam
      real(kind=rk)                         :: mx, result, epsabs, epsrel, abserr,ma, alphas
      integer(kind=ik)                      :: i, ier, neval, nd

      epsabs = 1e-7_rk
      epsrel = 1e-2_rk
      mx = params%mx
      ma = params%ma
      gam = 0.0_rk
      !if (sigma=="agff" .or. sigma=="afgf") then
        nd = size(params%alpha_s,2)
        call interp_linear(nd, params%alpha_s(1,:),params%alpha_s(2,:),T, alphas)
      !end if
      argsint%T=T

      select case(sigma)
      case("agff")
        do i=1, 9
          argsint%mf = mf(i)
          call qagi(kernel_agff,argsint,max(ma*ma,4.0_rk*mf(i)*mf(i)),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            !write(*,*) "Integral did not converge for the first tim ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
            call qags(kernel_agff,argsint,max(ma*ma,4.0_rk*mf(i)*mf(i)),&
                      1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + alpha_QED*qf(i)*qf(i)*ga*ggamma*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
        do i=4, 9
          argsint%mf = mf(i)
          call qagi(kernel_agff,argsint,max(1.0_rk,max(ma*ma,4.0_rk*mf(i)*mf(i))),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            call qags(kernel_agff,argsint,max(1.0_rk,max(ma*ma,4.0_rk*mf(i)*mf(i))),&
                      1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
            write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + 0.5_rk*alphas*ga*gg*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
      case("afgf")
        do i=1, 9
          argsint%mf = mf(i)
          call qagi(kernel_afgf,argsint,(ma+mf(i))*(ma+mf(i)),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            call qags(kernel_afgf,argsint,(ma+mf(i))*(ma+mf(i)),&
                      1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + alpha_QED*qf(i)*qf(i)*ga*gf(i)*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
        do i=4, 9
          argsint%mf = mf(i)
          call qagi(kernel_afgf,argsint,max(1.0_rk,(ma+mf(i))*(ma+mf(i))),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            call qags(kernel_afgf,argsint,max(1.0_rk,(ma+mf(i))*(ma+mf(i))),&
                      1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + 4.0_rk/3.0_rk*alphas*ga*gf(i)*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
      case("xxff")
        do i=1, 9
          argsint%mf = mf(i)
          argsint%nc = ncf(i)
          call qagi(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*mf(i)*mf(i)),&
                    1, epsabs, epsrel, result, abserr, neval, ier)
          if (ier > 0) then
            call qags(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*mf(i)*mf(i)),&
                      1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
          end if
          gam = gam + gDM*gDM*T/(32.0_rk*pi*pi*pi*pi)*result
        end do
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
      real(kind=rk)                         :: mx, result, epsabs, epsrel, abserr,ma,alphas
      integer(kind=ik)                      :: i, ier, neval, nd

      epsabs = 1e-20_rk
      epsrel = 1e-5_rk
      mx = params%mx
      ma = params%ma
      argsint%T=T
      if (sigma=="agff" .or. sigma=="afgf") then
        nd = size(params%alpha_s,2)
        call interp_linear(nd, params%alpha_s(1,:),params%alpha_s(2,:),T, alphas)
      end if
      select case(sigma)
        case("aaxx")
          if (T>0.5_rk) then
            call qagi(kernel_aaxx,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                  1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_aaxx,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                    1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
              if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
            end if
            sv = result/4.0_rk/(2.0_rk*T*ma*ma*ma*ma*bessK2(ma/T)* bessK2(ma/T))
          else
            call qagi(kernel_aaxx_series,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                  1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_aaxx_series,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                    1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
              if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
            end if
            sv = result
          end if
        case("xxaa")
          if (T>0.5_rk) then
            call qagi(kernel_xxaa,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                  1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_xxaa,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                    1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
              if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
            end if
            sv = result/4.0_rk/(2.0_rk*T*mx*mx*mx*mx* bessK2(mx/T)* bessK2(mx/T))
          else
            call qagi(kernel_xxaa_series,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                  1, epsabs, epsrel, result, abserr, neval, ier)
            if (ier > 0) then
              call qags(kernel_xxaa_series,argsint,max(4.0_rk*mx*mx,4.0_rk*ma*ma),&
                    1e10_rk, epsabs, epsrel, result, abserr, neval, ier)
              if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma
            end if
            sv = result
          end if
!        case("agff")
!          sv = 0.0_rk
!          do i=1, 9
!            argsint%mf = mf(i)
!            call qagi(kernel_agff,argsint,max(ma*ma,4.0_rk*mf(i)*mf(i)),&
!                1, epsabs, epsrel, result, abserr, neval, ier)
!            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
!            sv = sv + alpha_QED*qf(i)*qf(i)*result*1e10_rk/4.0_rk/(2.0_rk*T*ma*ma* bessK2(ma/T)*T*T*T)!????
!          end do
!          do i=4, 9
!            argsint%mf = mf(i)
!            call qagi(kernel_agff,argsint,max(1.0_rk,max(ma*ma,4.0_rk*mf(i)*mf(i))),&
!                1, epsabs, epsrel, result, abserr, neval, ier)
!            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
!            sv = sv + 0.5_rk*alphas*result*1e10_rk/4.0_rk/(2.0_rk*T*ma*ma* bessK2(ma/T)*T*T*T)!????
!          end do
!        case("afgf")
!          sv = 0.0_rk
!          do i=1, 9
!            argsint%mf = mf(i)
!            call qagi(kernel_afgf,argsint,(ma+mf(i))*(ma+mf(i)),&
!                1, epsabs, epsrel, result, abserr, neval, ier)
!            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
!            sv = sv + alpha_QED*qf(i)*qf(i)*result/4.0_rk/&
!                (2.0_rk*T*ma*ma* bessK2(ma/T)*mf(i)*mf(i)* bessK2(mf(i)/T))
!          end do
!          do i=4, 9
!            argsint%mf = mf(i)
!            call qagi(kernel_afgf,argsint,max(1.0_rk,(ma+mf(i))*(ma+mf(i))),&
!                1, epsabs, epsrel, result, abserr, neval, ier)
!            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
!            sv = sv +  4.0_rk/3.0_rk*alphas*result/4.0_rk/(2.0_rk*T*ma*ma* &
!                bessK2(ma/T)*mf(i)*mf(i)* bessK2(mf(i)/T))
!          end do
!        case("xxff")
!          sv = 0.0_rk
!          do i=1,9
!            argsint%mf = mf(i)
!            argsint%nc = ncf(i)
!            call qagi(kernel_xxff,argsint,max(4.0_rk*mx*mx,4.0_rk*mf(i)*mf(i)),&
!                      1, epsabs, epsrel, result, abserr, neval, ier)
!            if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "sigma=", sigma, " fermion=", i
!            sv = sv + result/4.0_rk/(2.0_rk*T*mx*mx*mf(i)*mf(i)* &
!                      bessK2(mx/T)* bessK2(mf(i)/T))
!          end do
        case default
          write(*,*) "Error! x section", sigma, "not implemented"
      end select
    end subroutine sigmav

    real(kind=rk) function kernel_agff( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint
        kernel_agff = sigma_agff(s,argsint%mf,argsint%ma,argsint%g)*4.0_rk*&
                      F(s,argsint%ma,0.0_rk)*F(s,argsint%ma,0.0_rk)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel_agff
    real(kind=rk) function kernel_ffag( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint
        kernel_ffag = sigma_ffag(s,argsint%mf,argsint%ma,argsint%g)*4.0_rk*&
                      F(s,argsint%mf,argsint%mf)*F(s,argsint%mf,argsint%mf)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel_ffag
    real(kind=rk) function kernel_xxff( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint

        kernel_xxff = sigma_xxff(s,argsint%mf,argsint%nc,argsint%mx,argsint%ma,argsint%g)*4.0_rk*&
                      F(s,argsint%mx,argsint%mx)*F(s,argsint%mx,argsint%mx)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel_xxff
    real(kind=rk) function kernel_afgf( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint
        kernel_afgf = sigma_afgf(s,argsint%mf,argsint%ma,argsint%g)*4.0_rk*&
                      F(s,argsint%ma,argsint%mf)*F(s,argsint%ma,argsint%mf)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel_afgf
    real(kind=rk) function kernel_aaxx( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint
        kernel_aaxx = sigma_aaxx(s,argsint%mx,argsint%ma,argsint%g)*4.0_rk*&
                      F(s,argsint%ma,argsint%ma)*F(s,argsint%ma,argsint%ma)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
        return
    end function kernel_aaxx

    real(kind=rk) function kernel_xxaa( s, argsint )
      implicit none
        real(kind=rk), intent(in)           :: s
        type (type_argsint), intent(in)     :: argsint

        kernel_xxaa = sigma_xxaa(s,argsint%mx,argsint%ma,argsint%g)*4.0_rk*&
                      F(s,argsint%mx,argsint%mx)*F(s,argsint%mx,argsint%mx)/sqrt(s) * bessK1( sqrt(s)/argsint%T)

    end function kernel_xxaa

    real(kind=rk) function F(s,m1,m2)
      implicit none
        real(kind=rk), intent(in)  :: s, m1, m2
        F = 0.5_rk*sqrt((s-(m1+m2)*(m1+m2))*(s-(m1-m2)*(m1-m2)))
        return
    end function F

    real(kind=rk) function kernel_xxaa_series( s, argsint )
      implicit none
      real(kind=rk), intent(in)           :: s
      type (type_argsint), intent(in)     :: argsint
      real(kind=rk)                       :: mx, T

      mx = argsint%mx
      T = argsint%T

      kernel_xxaa_series = sigma_xxaa(s,mx,argsint%ma,argsint%g)*&
                    F(s,mx,mx)*F(s,mx,mx)/sqrt(s)*&
                    (exp((2.0_rk*mx - sqrt(s))/T)* mx* &
                    (1140.0_rk*s*T*T -60.0_rk*mx*T*(8.0_rk*s + 3.0_rk*sqrt(s)*T)&
                    +mx*mx*(128.0_rk*s + 48.0_rk*sqrt(s)*T - 15.0_rk*T*T)))&
                    /(128.0_rk*mx*mx*mx*mx*mx*sqrt(2.0_rk*pi)*s**1.25_rk*T**1.5_rk)

    end function kernel_xxaa_series

    real(kind=rk) function kernel_aaxx_series( s, argsint )
      implicit none
      real(kind=rk), intent(in)           :: s
      type (type_argsint), intent(in)     :: argsint
      real(kind=rk)                       :: ma, T

      ma = argsint%ma
      T = argsint%T

      kernel_aaxx_series = sigma_aaxx(s,argsint%mx,ma,argsint%g)*&
                    F(s,ma,ma)*F(s,ma,ma)/sqrt(s)*&
                    (exp((2.0_rk*ma - sqrt(s))/T)* ma* &
                    (1140.0_rk*s*T*T -60.0_rk*ma*T*(8.0_rk*s + 3.0_rk*sqrt(s)*T)&
                    +ma*ma*(128.0_rk*s + 48.0_rk*sqrt(s)*T - 15.0_rk*T*T)))&
                    /(128.0_rk*ma*ma*ma*ma*ma*sqrt(2.0_rk*pi)*s**1.25_rk*T**1.5_rk)

    end function kernel_aaxx_series

    !include "rhs_boltzmann.f90"
    !include "rhs_region3a2.f90"
    !include "region3a_log.f90"
    include "rhs_contributions.f90"
    include "region3aeq.f90"
    include "region3a_in_n.f90"
    include "RK4.f90"
end module
