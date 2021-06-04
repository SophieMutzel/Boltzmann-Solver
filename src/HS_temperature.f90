
    real(kind=rk) function Taroot(Ta,T,params,rhoprime)
      implicit none
      real(kind=rk), intent(in)           :: Ta, T,rhoprime
      type (type_params), intent(in)      :: params
      integer(kind=ik)                    :: nd, nr
      real(kind=rk)                       :: gHS, gSM, rar, rhoa, rhoDM

      !nd = size(params%geff_HS,2)
      ! DOFS
      !call test_bezier(params%geff_HS(1,:),params%geff_HS(2,:),log10(Ta),gHS,nd-1,params%A,params%B)
      !call geffSM(T,params,gSM)
      gSM = geff_rho(T)
      ! Axion temperature
      !nr = size(params%rhoa_rho,2)
      !call interp_linear(nr, params%rhoa_rho(1,:),params%rhoa_rho(2,:),T, rar)
      rhoa = rhoeq(Ta,params%ma,ga)
      rhoDM = rhoeq(Ta,params%mx,gDM)
      !Taroot = sqrt(params%gaff(1)*sqrt(( 10.75_rk/(1.0_rk) *rar )))*T !- Ta
      Taroot = (rhoa+rhoDM)/(gSM*T*T*T*T*pi*pi/30.0_rk)-params%gaff(1)*params%gaff(1)*rhoprime
      !if (T>params%mx) then
      !  Taroot = rhoa+rhoDM - rhoprime
      !else
      !  Taroot = rhoa+rhoDM*facrhoDM -rhoprime
      !end if
    end function Taroot

!    real(kind=rk) function Tanew(T,params,rhoprime,facrhoDM)
!      implicit none
!      real(kind=rk), intent(in)               :: T
!      type (type_params), intent(in)          :: params
!      real(kind=rk), dimension(:), intent(in) :: facrhoDM,rhoprime
!      integer(kind=ik)                        :: i
!
!      do i=1,params%N
!        Tanew=rtbis(Taroot,1e-5_rk,100000.0_rk,1e-10_rk,params,T,rhoprime(i),facrhoDM(i))
!      end do
!    end function Tanew
!
    real(kind=rk) function Ta(T,params,rhoa)
      implicit none
      real(kind=rk), intent(in)           :: T, rhoa
      type (type_params), intent(in)      :: params
      Ta = rtbis(Taroot,1e-4_rk,10000.0_rk,1e-10_rk,params,T,rhoa)
      !Ta = Taroot(T,T,params)!
    end function Ta

    FUNCTION rtbis(func,x1,x2,xacc,params, T,rhoprime)
    IMPLICIT NONE
    REAL(kind=rk), INTENT(IN) :: x1,x2,xacc, T,rhoprime
    type (type_params), intent(in)      :: params
    REAL(kind=rk) :: rtbis
    real ( kind = rk ), external :: func
    INTEGER(kind=ik), PARAMETER :: MAXIT=80
    !Using bisection, find the root of a function func known to lie between x1 and x2. The root, returned as rtbis, will be refined until its accuracy is ±xacc.
    !Parameter: MAXIT is the maximum allowed number of bisections.
    INTEGER(kind=ik) :: j
    REAL(kind=rk) :: dx,f,fmid,xmid
    fmid=func(x2,T,params,rhoprime)
    f=func(x1,T,params,rhoprime)
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
      fmid=func(xmid,T,params,rhoprime)
      if (fmid <= 0.0) rtbis=xmid
      if (abs(dx) < xacc .or. fmid == 0.0) RETURN
    end do
    call abort_it("rtbis: too many bisections")
    END FUNCTION rtbis
