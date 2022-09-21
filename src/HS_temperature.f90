! routines needed to calculate the HS temperature for intial condition
    real(kind=rk) function Taroot(Ta,T,params,rhoprime)
      implicit none
      real(kind=rk), intent(in)           :: Ta, T,rhoprime
      type (type_params), intent(in)      :: params
      real(kind=rk)                       :: rar, rhoa, rhoDM

      rhoa = rhoeq(Ta,params%ma,ga)
      rhoDM = rhoeq(Ta,params%mx,gDM)
      Taroot = (rhoa+rhoDM)/rho_SM(T)-params%gaff(1)*params%gaff(1)*rhoprime
    end function Taroot

    real(kind=rk) function Ta(T,params,rhoa)
      implicit none
      real(kind=rk), intent(in)           :: T, rhoa
      type (type_params), intent(in)      :: params
      Ta = rtbis(Taroot,1e-5_rk,1e5_rk,1e-5_rk,params,T,rhoa)
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

    ! we don't use these anymore but assume T_ALPs=T_SM in sequential freeze-in
    real(kind=rk) function Ta_seq_fi(T,params,rhooverna)
      implicit none
      real(kind=rk), intent(in)           :: T, rhooverna
      type (type_params), intent(in)      :: params
      Ta_seq_fi = rtbis(Taroot_seq_fi,1e-6_rk,100.0_rk,1e-10_rk,params,T,rhooverna)
    end function Ta_seq_fi

    real(kind=rk) function Taroot_seq_fi(Ta,T,params,rhoovern)
      implicit none
      real(kind=rk), intent(in)           :: Ta, T,rhoovern
      type (type_params), intent(in)      :: params

      Taroot_seq_fi = rhoeqneq(Ta,params%ma)-rhoovern
    end function Taroot_seq_fi
