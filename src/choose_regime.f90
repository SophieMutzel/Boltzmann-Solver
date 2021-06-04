subroutine choose_regime(params, argsint)
  implicit none
  type (type_params), intent(inout)  :: params
  type (type_argsint), intent(inout) :: argsint
  real(kind=rk)                      :: mx, ma, abserr, rar, Tmx
  real(kind=rk)                      :: neqazp, sv_aaxx, gam_agff, ffa, gam_xxff
  real(kind=rk)                      :: SM_axion, SM_DM, axion_DM
  real(kind=rk), dimension(2,5)      :: Tprime
  integer(kind=ik)                   :: i, ier, neval

  Tprime(2,:) = (/0.1_rk,1.0_rk,10.0_rk,100.0_rk,500.0_rk/)
  mx = params%mx
  ma = params%ma
  ! SM axion interaction
  call gamma_r_new( ma, params, argsint, "agffth", gam_agff )
  gam_agff = gam_agff*params%gaff(1)*params%gaff(1)
  ! inverse decay ff->a
  ffa = gammav(ma, params, "affth")*params%gaff(1)*params%gaff(1)
  !call gamma_r_new( T, params, argsint, "afgf", gam_afgf(1) )
  ! SM DM interaction
  call gamma_r_new( mx, params, argsint, "xxffth", gam_xxff )
  gam_xxff = gam_xxff*params%gaxx(1)*params%gaff(1)*params%gaxx(1)*params%gaff(1)

  ! check equilibrium SM<->axions
  SM_axion = (gam_agff+ffa)/neq(ma,ma,ga)/Hub(ma,rhoeq(ma,mx,gDM)+rhoeq(ma,ma,ga))
  ! check equilibrium SM<->DM
  SM_DM = gam_xxff/neq(mx,mx,gDM)/Hub(mx,rhoeq(mx,mx,gDM)+rhoeq(mx,ma,ga))

  if (( SM_axion > 1.0_rk ) .and. ( SM_DM > 1.0_rk )) then
      params%regime = "freeze-out"
  else
    if ( SM_axion > 1.0_rk ) then
      params%regime = "freeze-in"
    else
      allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
      argsint%drhoa = params%drhoa
      do i=1,5
        call qags(rhop_over_rho,argsint,T_RH,&
                  Tprime(2,i), 1e-5_rk, 1e-5_rk, rar, abserr, neval, ier)
        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for T=", Tprime(2,i)
        Tprime(1,i) = Ta(Tprime(2,i),params,rar)
      end do
      deallocate(argsint%drhoa)
      ! find T for which T' is mx
      call interp_linear(5, Tprime(1,:),Tprime(2,:),mx, Tmx)
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
    end if
  end if

end subroutine choose_regime
