subroutine choose_regime(params, argsint)
  implicit none
  type (type_params), intent(inout)  :: params
  type (type_argsint), intent(inout) :: argsint
  real(kind=rk)                      :: mx, ma, abserr, rar, Tmx
  real(kind=rk)                      :: neqazp, sv_aaxx, gam_agff, ffa, gam_xxff,gam_afgf
  real(kind=rk)                      :: SM_axion, SM_DM, axion_DM
  integer(kind=ik)                   :: i, ier, neval, nT=20
  real(kind=rk), allocatable         :: Tprime(:,:)

  mx = params%mx
  ma = params%ma
  ! SM axion interaction
  call gamma_r_new( mx, argsint, "agffth", gam_agff )
  gam_agff = gam_agff*params%gaff(1)*params%gaff(1)
  ! inverse decay a->ff
  ffa = gammav(mx, argsint, "affth")*params%gaff(1)*params%gaff(1)

  call gamma_r_new( mx, argsint, "afgfth", gam_afgf )
  gam_afgf = gam_afgf*params%gaff(1)*params%gaff(1)
  ! SM DM interaction
  call gamma_r_new( mx, argsint, "xxffth", gam_xxff )
  gam_xxff = gam_xxff*params%gaxx(1)*params%gaff(1)*params%gaxx(1)*params%gaff(1)

  ! check equilibrium SM<->axions
  SM_axion = (gam_agff/neq(mx,ma,ga)+2.0_rk*gam_afgf/neq(mx,ma,ga)+ffa)/Hub(mx,rhoeq(mx,mx,gDM)+rhoeq(mx,ma,ga))
  ! check equilibrium SM<->DM
  SM_DM = gam_xxff/neq(mx,mx,gDM)/Hub(mx,rhoeq(mx,mx,gDM)+rhoeq(mx,ma,ga))

  if (( SM_axion > 1.0_rk ) .and. ( SM_DM > 1.0_rk )) then
      params%regime = "freeze-out"
  else
    if ( SM_axion > 1.0_rk ) then
      params%regime = "freeze-in"
    else
      allocate(Tprime(2,nT))
      !Tprime(2,:) = (/0.1_rk,1.0_rk,10.0_rk,100.0_rk,500.0_rk/)
      call linspace(-2.0_rk, log10(T_RH)-0.2_rk, Tprime(2,:))
      Tprime(2,:) = 10.0_rk**Tprime(2,:)
      allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
      argsint%drhoa = params%drhoa
      do i=1,nT
        call qags(rhop_over_rho,argsint,T_RH,&
                  Tprime(2,i), 1e-5_rk, 1e-5_rk, rar, abserr, neval, ier)
        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for T=", Tprime(2,i)
        Tprime(1,i) = Ta(Tprime(2,i),params,rar)
      end do
      deallocate(argsint%drhoa)
      ! find T for which T' is mx
      call interp_linear(nT, Tprime(1,:),Tprime(2,:),mx, Tmx)
      ! neq,a(z')
      neqazp = neq(mx, ma, ga)
      ! HS interaction
      call sigmav( mx, params, argsint, "aaxx", sv_aaxx )
      sv_aaxx = sv_aaxx*params%gaxx(1)*params%gaxx(1)*params%gaxx(1)*params%gaxx(1)
      ! check equilibrium DM<->axions
      axion_DM = sv_aaxx*neqazp/Hub(Tmx,rhoeq(mx,mx,gDM)+rhoeq(mx,ma,ga))
      if ( axion_DM > 1.0_rk) then
        params%regime = "reannihilation"
      else
        params%regime = "seq-freeze-in"
      end if
      deallocate(Tprime)
    end if
  end if
end subroutine choose_regime
