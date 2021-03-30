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

        Hub = 1.67 * geff(T,params) * T* T / Mpl
        return
    end function Hub

    real(kind=rk) function ent( T, params )
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

    real(kind=rk) function Taroot(Ta,T,params)
      implicit none
      real(kind=rk), intent(in)           :: Ta, T
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
      rhoa = ga*params%ma*params%ma*Ta*(params%ma*bessK1(params%ma/Ta)+3.0_rk*Ta*bessK2(params%ma/Ta))
      rhoDM = gDM*params%mx*params%mx*Ta*(params%mx*bessK1(params%mx/Ta)+3.0_rk*Ta*bessK2(params%mx/Ta))
      !Taroot = sqrt(params%gaff(1)*sqrt(( 10.75_rk/(1.0_rk) *rar )))*T !- Ta
      Taroot = 15.0_rk*(rhoa+rhoDM)/(gSM*T*T*T*T*pi*pi*pi*pi)-params%gaff(1)*params%gaff(1)*rar

    end function Taroot

    real(kind=rk) function Ta(T,params)
      implicit none
      real(kind=rk), intent(in)           :: T
      type (type_params), intent(in)      :: params
      Ta = rtbis(Taroot,1e-3_rk,10000.0_rk,1e-10_rk,params,T)
      !Ta = Taroot(T,T,params)!
    end function Ta

    FUNCTION rtbis(func,x1,x2,xacc,params, T)
    IMPLICIT NONE
    REAL(kind=rk), INTENT(IN) :: x1,x2,xacc, T
    type (type_params), intent(in)      :: params
    REAL(kind=rk) :: rtbis
    real ( kind = rk ), external :: func
    INTEGER(kind=ik), PARAMETER :: MAXIT=80
    !Using bisection, find the root of a function func known to lie between x1 and x2. The root, returned as rtbis, will be refined until its accuracy is ±xacc.
    !Parameter: MAXIT is the maximum allowed number of bisections.
    INTEGER(kind=ik) :: j
    REAL(kind=rk) :: dx,f,fmid,xmid
    fmid=func(x2,T,params)
    f=func(x1,T,params)
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
      fmid=func(xmid,T,params)
      if (fmid <= 0.0) rtbis=xmid
      if (abs(dx) < xacc .or. fmid == 0.0) RETURN
    end do
    call abort_it("rtbis: too many bisections")
    END FUNCTION rtbis

    include "initial_conditions.f90"
end module
