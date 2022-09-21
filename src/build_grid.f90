! build grid for phase diagram
subroutine build_grid(params, argsint, gaff_tot, gaxx_tot)
implicit none
type (type_params), intent(inout)       :: params
type (type_argsint), intent(inout)      :: argsint
real(kind=rk), intent(out), allocatable :: gaff_tot(:), gaxx_tot(:)
real(kind=rk), allocatable              :: gaxx(:), gaff(:), gaff_vec(:), gaxx_vec(:), ma_QCD(:),gaff_QCD(:),neqx(:)!, Tprime_vec(:)
real(kind=rk)                           :: T_start, abserr, rar, Tmx, gaff_low, gaff_up, ma_low
real(kind=rk)                           :: aaxx, epsabs, epsrel
integer(kind=ik)                        :: i, ier, neval, nT=21, n_gaxx=10, n_gaff, j, k, pos
real(kind=rk)                           :: Tprime, delta=0.5_rk, gaff_up_QCD,neqobs


! lower bound gaff
  allocate(argsint%drhoa(2,size(params%drhoa(1,:))),gaxx(nT),gaff(nT),neqx(nT))
  argsint%g = 1.0_rk
  epsabs=1e-30_rk
  epsrel=1e-30_rk
  argsint%drhoa = params%drhoa
  call linspace(-8.0_rk, -15.0_rk, gaff)

  !allocate(ma(nT))!,Tprime_vec(nT))
!  call linspace(-3.0_rk,2.0_rk,ma)
!  T_start=QCDcut
!  do i=1,nT
!    do k=1,nT
!    !Tprime=10000.0_rk
!    !ma=-0.1_rk
!    !params%mx=1.0_rk
!    !do while (Tprime>params%mx/20.0_rk .and. ma>-3.0_rk)
!      params%gaff(1)=10.0_rk**gaff(i)
!      params%ma=10.0**ma(k)
!      params%mx=params%ma*10.0_rk
!      call qags(rhop_over_rho,argsint,T_RH,&
!            T_start, epsabs, epsrel, rar, abserr, neval, ier)
!      !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
!      Tprime = Ta(T_start,params,rar)
!      !ma = ma-0.005_rk
!      write(*,*) params%gaff(1), params%ma, Tprime
!    end do
!    !call interp_linear(nT, ma,Tprime_vec,params%mx, ma_low)
!  end do
!
!  stop
  allocate(ma_QCD(nT),gaff_QCD(nT))

  ma_QCD=(/-2.50559, -2.32372, -2.14772, -1.97642, -1.79896, -1.62184, &
          -1.44856, -1.27535, -1.09694, -0.922434, -0.752097, -0.572615, &
          -0.396722, -0.22523, -0.0487343, 0.128436, 0.301483, 0.474775, &
          0.653136, 0.827641, 0.997946/)
  gaff_QCD=(/-15., -14.65, -14.3, -13.95, -13.6, -13.25, -12.9, -12.55, -12.2, &
          -11.85, -11.5, -11.15, -10.8, -10.45, -10.1, -9.75, -9.4, -9.05, &
          -8.7, -8.35, -8./)
  do i=1,nT
    params%gaff(1)=10.0_rk**gaff(i)
    params%z_start = log10(params%ma)-1.5_rk
    Tprime=1000000.0_rk
    do while (Tprime>2.0_rk*params%mx)
      params%z_start = params%z_start+0.2_rk
      T_start = params%mx/10**params%z_start
      call qags(rhop_over_rho,argsint,T_RH,&
            T_start, epsabs, epsrel, rar, abserr, neval, ier)
      if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
      Tprime = Ta(T_start,params,rar)
    end do
  !  call sigmav( Tprime, params, argsint, "aaxx", aaxx )
  !  gaxx(i) = log10(sqrt(sqrt(Hub(T_start,rhoeq(Tprime,params%mx,gDM)+&
  !        rhoeq(Tprime,params%ma,ga))/aaxx/neq(Tprime,params%ma,ga))))
    if (Tprime>T_start) then
      gaff_up = gaff(i)
      do while (Tprime>T_start)
        gaff_up = gaff_up - 0.1_rk
        gaff(i) = gaff_up
        params%gaff(1) = 10.0_rk**gaff_up
        params%z_start = log10(params%ma)-1.5_rk
        Tprime=1000000.0_rk
        do while (Tprime>params%mx)
          params%z_start = params%z_start+0.2_rk
          T_start = params%mx/10**params%z_start
          call qags(rhop_over_rho,argsint,T_RH,&
                T_start, epsabs, epsrel, rar, abserr, neval, ier)
          if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
          Tprime = Ta(T_start,params,rar)
        end do
      end do
    end if
  end do
!  if ( maxval(gaxx)*minval(gaxx) < 0.0_rk ) then
!!    do i=1,size(gaxx)
!!      if (gaxx(i)>0.0_rk) then
!!        pos = i
!!        exit
!!      end if
!!    end do
!    call interp_linear(nT, gaxx,gaff,0.0_rk, gaff_low)
!  else
!    if ( maxval(gaxx) < 0.0_rk ) then
!      gaff_low = -15.0_rk
!    else
!      call abort_it("gaxx larger than 1 for these masses! ma= "//float2str(params%ma)//" mx= "//float2str(params%mx))
!    end if
!  end if

  do i=1,nT
    params%gaff(1)=10.0_rk**gaff(i)
    params%z_start = log10(params%ma)-1.5_rk
    Tprime=1000000.0_rk
    do while (Tprime>params%mx)
      params%z_start = params%z_start+0.1_rk
      T_start = params%mx/10**params%z_start
      call qags(rhop_over_rho,argsint,T_RH,&
            T_start, epsabs, epsrel, rar, abserr, neval, ier)
      !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
      Tprime = Ta(T_start,params,rar)
    end do
    neqx(i) = neq(Tprime,params%mx,gDM)
    neqobs = omegax*rhocrit/params%mx/2.0_rk/s0*ent(T_start,params)
    if (neqx(i)>neqobs) then
      gaff_low = gaff(i)
    end if
  end do
  !call interp_linear(nT, neqx,gaff,neqobs, gaff_low)
  gaff_low = gaff_low + 0.1_rk!-13.5
  call interp_linear(nT, ma_QCD,gaff_QCD,log10(params%ma), gaff_up_QCD)
  gaff_up = min(maxval(gaff),gaff_up_QCD)
  !gaff_up = maxval(gaff)
  deallocate(gaxx,gaff)
  n_gaff = int((gaff_up-gaff_low)/delta)+1
  allocate(gaxx_tot(n_gaxx*n_gaff),gaff_tot(n_gaxx*n_gaff),gaff_vec(n_gaff),gaxx_vec(n_gaxx))
  call linspace(gaff_low, gaff_up, gaff_vec)
  gaff_vec=10.0_rk**gaff_vec
  call linspace(-3.5_rk, 0.0_rk, gaxx_vec)
  gaxx_vec=10.0_rk**gaxx_vec
  k = 0
  do i=1,n_gaff
    do j=1,n_gaxx
      k = k + 1
      gaxx_tot(k) = gaxx_vec(j)
      gaff_tot(k) = gaff_vec(i)
    end do
  end do
  deallocate(gaxx_vec,gaff_vec,argsint%drhoa)
  !gaff_up = sqrt(1.0_rk/rar)
!  allocate(Tprime(2,nT))
!  call linspace(-1.0_rk, 2.7_rk, Tprime(2,:))
!  Tprime(2,:) = 10.0_rk**Tprime(2,:)
!  allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
!  argsint%drhoa = params%drhoa
!  do i=1,nT
!    call qags(rhop_over_rho,argsint,T_RH,&
!              Tprime(2,i), 1e-5_rk, 1e-5_rk, rar, abserr, neval, ier)
!    if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for T=", Tprime(2,i)
!    Tprime(1,i) = Ta(Tprime(2,i),params,rar)
!  end do
!  deallocate(argsint%drhoa)
!  ! find T for which T' is mx
!  call interp_linear(nT, Tprime(1,:),Tprime(2,:),mx, Tmx)
!  ! neq,a(z')
!  neqazp = neq(mx, ma, ga)
!  ! HS interaction
!  call sigmav( mx, params, argsint, "aaxx", sv_aaxx )
!  sv_aaxx = sv_aaxx*params%gaxx(1)*params%gaxx(1)*params%gaxx(1)*params%gaxx(1)
!  ! check equilibrium DM<->axions
!  axion_DM = sv_aaxx*neqazp/Hub(Tmx,rhoeq(mx,mx,gDM)+rhoeq(mx,ma,ga))

!  argsint%g = 1.0_rk
!  allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
!  argsint%drhoa = params%drhoa
!  do j=-3,3,0.1
!    ma=10.0_rk**j
!    argsint%ma=ma
!    argsint%mx=10.0_rk*argsint%ma
!    epsabs=1e-30_rk
!    epsrel=1e-30_rk
!    call qags(rhop_over_rho,argsint,T_RH,&
!          T, epsabs, epsrel, rar, abserr, neval, ier)
!    write(*,*) j, sqrt(1.0_rk/rar)
!    !sqrt(0.5_rk/rar/rho_SM(0.001_rk)*7.0_rk/8.0_rk*ggamma*pi*pi/30.0_rk*1e-12_rk)
!  end do
!  stop
end subroutine build_grid
