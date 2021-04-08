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

        Hub = 1.67_rk * geff(T,params) * T* T / Mpl
        return
    end function Hub

    real(kind=rk) function ent( T, params )
      implicit none
        real(kind=rk), intent(in)           :: T
        type (type_params), intent(in)      :: params

        ent = 2.0_rk * pi* pi /45.0_rk * heff(T,params) * T*T*T
    end function ent

    real(kind=rk) function neq( T, m, g )
      implicit none
      real(kind=rk), intent(in)           :: T, m, g

      neq = 0.5_rk * g * T * m*m * bessK2(m/T) /pi/pi

    end function neq

    real(kind=rk) function rhoeq( T, m, g )
      implicit none
      real(kind=rk), intent(in)           :: T, m, g

      rhoeq = 0.5_rk * g *m*m*T*(m*bessK1(m/T)+3.0_rk*T*bessK2(m/T)) /pi/pi
    end function rhoeq

    real(kind=rk) function peq( T, m, g )
      implicit none
      real(kind=rk), intent(in)           :: T, m, g

      peq = 0.5_rk * g *m*m*T*T*bessK2(m/T)/pi/pi

    end function peq

    real(kind=rk) function drhoeq( T, m, g )
      ! d(rho_eq(T)/dT
      implicit none
      real(kind=rk), intent(in)           :: T, m, g
      drhoeq = g*m/(2.0_rk*pi*pi)*((m*m*m/T + 12.0_rk*m*T)*bessK0(m/T) + (5.0_rk*m*m + 24.0_rk*T*T)*BessK1(m/T))
    end function drhoeq

    real(kind=rk) function drhoeqneq( T, m, g )
      ! d(rho_eq(T)/neq(T))/dT
      implicit none
      real(kind=rk), intent(in)           :: T, m, g
      real(kind=rk)                       :: bk1, bk2

      !bk1 = bessK1(m/T)
      !bk2 = bessK2(m/T)
      !drhoeqneq = -((m*m* bk1*bk1 - bk2*(m*m*bessK0(m/T) + (m*m + 6.0_rk* T*T)*bk2) + &
      !            m*m*bk1*bessk_s(3,m/T))/(2.0_rk* T*T* bk2*bk2))
      drhoeqneq = 3.0_rk/2.0_rk + (15.0_rk*T)/(4.0_rk*m) - (45.0_rk*T*T)/(8.0_rk*m*m)
    end function drhoeqneq

    real(kind=rk) function Taroot(Ta,T,params,rhoprime,facrhoDM)
      implicit none
      real(kind=rk), intent(in)           :: Ta, T,rhoprime,facrhoDM
      type (type_params), intent(in)      :: params
      integer(kind=ik)                    :: nd, nr
      real(kind=rk)                       :: gHS, gSM, rar, rhoa, rhoDM

      !nd = size(params%geff_HS,2)
      ! DOFS
      !call test_bezier(params%geff_HS(1,:),params%geff_HS(2,:),log10(Ta),gHS,nd-1,params%A,params%B)
      !call geffSM(T,params,gSM)
      gSM = geff_rho(T)
      ! Axion temperature
      nr = size(params%rhoa_rho,2)
      call interp_linear(nr, params%rhoa_rho(1,:),params%rhoa_rho(2,:),T, rar)
      rhoa = rhoeq(Ta,params%ma,ga)
      rhoDM = rhoeq(Ta,params%mx,gDM)
      !Taroot = sqrt(params%gaff(1)*sqrt(( 10.75_rk/(1.0_rk) *rar )))*T !- Ta
      Taroot = (rhoa+rhoDM)/(gSM*T*T*T*T*pi*pi/30.0_rk)-params%gaff(1)*params%gaff(1)*rar
      !if (T>params%mx) then
      !  Taroot = rhoa+rhoDM - rhoprime
      !else
      !  Taroot = rhoa+rhoDM*facrhoDM -rhoprime
      !end if
    end function Taroot

    real(kind=rk) function Tanew(T,params,rhoprime,facrhoDM)
      implicit none
      real(kind=rk), intent(in)               :: T
      type (type_params), intent(in)          :: params
      real(kind=rk), dimension(:), intent(in) :: facrhoDM,rhoprime
      integer(kind=ik)                        :: i

      do i=1,params%N
        Tanew=rtbis(Taroot,1e-5_rk,100000.0_rk,1e-10_rk,params,T,rhoprime(i),facrhoDM(i))
      end do
    end function Tanew

    real(kind=rk) function Ta(T,params)
      implicit none
      real(kind=rk), intent(in)           :: T
      type (type_params), intent(in)      :: params
      Ta = rtbis(Taroot,1e-3_rk,1000.0_rk,1e-10_rk,params,T,1.0_rk,1.0_rk)
      !Ta = Taroot(T,T,params)!
    end function Ta

    FUNCTION rtbis(func,x1,x2,xacc,params, T,rhoprime,facrhoDM)
    IMPLICIT NONE
    REAL(kind=rk), INTENT(IN) :: x1,x2,xacc, T,facrhoDM,rhoprime
    type (type_params), intent(in)      :: params
    REAL(kind=rk) :: rtbis
    real ( kind = rk ), external :: func
    INTEGER(kind=ik), PARAMETER :: MAXIT=80
    !Using bisection, find the root of a function func known to lie between x1 and x2. The root, returned as rtbis, will be refined until its accuracy is ±xacc.
    !Parameter: MAXIT is the maximum allowed number of bisections.
    INTEGER(kind=ik) :: j
    REAL(kind=rk) :: dx,f,fmid,xmid
    fmid=func(x2,T,params,rhoprime,facrhoDM)
    f=func(x1,T,params,rhoprime,facrhoDM)
    if (f*fmid >= 0.0) call abort_it("rtbis: root must be bracketed")
    ! Orient the search so that f>0 lies at x+dx.
    if (f < 0.0) then
      rtbis=x1
      dx=x2-x1
    else
      rtbis=x2
      dx=x1-x2
    end if

    do j=1,MAXIT
      dx=dx*0.5_rk
      !Bisection loop.
      xmid=rtbis+dx
      fmid=func(xmid,T,params,rhoprime,facrhoDM)
      if (fmid <= 0.0) rtbis=xmid
      if (abs(dx) < xacc .or. fmid == 0.0) RETURN
    end do
    call abort_it("rtbis: too many bisections")
    END FUNCTION rtbis

    include "initial_conditions.f90"
end module
