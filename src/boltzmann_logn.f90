! rhs of Boltzmann equations, solving log10(n)
subroutine boltzmann_logn( N, lz, Y, Ynew, params, argsint )
  implicit none
  real(kind=rk), intent(in)                   :: lz
  real(kind=rk), dimension(N), intent(in)     :: Y
  real(kind=rk), dimension(N), intent(inout)  :: Ynew
  type (type_params), intent(in)              :: params
  type (type_argsint), intent(inout)          :: argsint
  integer(kind=ik), intent(in)                :: N
  real(kind=rk), dimension(nrhs,params%N)     :: q, rhs, q10
  real(kind=rk)                               :: mx, ma, T, mf, nc, s, H, drhoa
  real(kind=rk)                               :: result
  integer(kind=ik)                            :: ier, i, neval,nd,nz
  real(kind=rk), dimension(params%N)          :: Tp, neqzp,neqazp, neq10, neqa10
  real(kind=rk), dimension(params%N)          :: gam_agff, gam_afgf, gam_xxff, ffa,sv_aaxx, sv_xxaa
  real(kind=rk)                               :: rhoeqneqa,rhoeqneqDM
  real(kind=rk)                               :: rhoeqaTp, rhoeqDMTp, peqaTp, peqDMTp,rhoplusp

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  ma = params%ma
  T = mx/10**lz
  q10 = 10.0_rk**q
  if (params%regime=="reannihilation") then
    Tp(:) = q(3,:)
    nd = size(params%drhoa,2)
    call interp_linear(nd, params%drhoa(1,:),params%drhoa(2,:),log10(T), drhoa)
    ! inverse decays ff->a
    drhoa = drhoa + drho_decay(T, params%ma, "ffath")
  else
    Tp = T
  end if
  ! rho,eq,a(T')
  rhoeqaTp = rhoeq(Tp(1),ma,ga)
  ! rhoeq,DM(T')
  rhoeqDMTp = rhoeq(Tp(1),mx,gDM)
  ! Hubble function
  H = Hub( T, rhoeqaTp+rhoeqDMTp )
  do i=1,params%N
    ! neq,DM(z')
    neqzp(i) = neq(Tp(i), mx, gDM)
    ! neq,a(z')
    neqazp(i) = neq(Tp(i), ma, ga)
    ! HS interaction
    call sigmav( Tp(i), params, argsint, "aaxx", sv_aaxx(i) )
    call sigmav( Tp(i), params, argsint, "xxaa", sv_xxaa(i) )
    sv_aaxx(i) = sv_aaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    sv_xxaa(i) = sv_xxaa(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    ! SM axion interaction
    call gamma_r_new( T, argsint, "agffth", gam_agff(i) )
    gam_agff(i) = gam_agff(i)*params%gaff(i)*params%gaff(i)

    nz = size(params%gam_afgf,2)
    call interp_linear(nz, params%gam_afgf(1,:),params%gam_afgf(2,:),lz, gam_afgf(i))
    !call gamma_r_new( T, argsint, "afgfth", gam_afgf(i) )
    gam_afgf(i) = gam_afgf(i)*params%gaff(i)*params%gaff(i)
    ! SM DM interaction
    call gamma_r_new( T, argsint, "xxffth", gam_xxff(i) )
    gam_xxff(i) = gam_xxff(i)*params%gaxx(i)*params%gaff(i)*params%gaxx(i)*params%gaff(i)
    ! inverse decay ff->a
    ffa(i) = gammav(T, argsint, "affth")*params%gaff(i)*params%gaff(i)
  end do

  select case (params%regime)
  case ("reannihilation")
    if (Tp(1)>ma) then
      rhs(1,:) =  - 3.0_rk +  ((-sv_xxaa*q10(1,:)+sv_aaxx* 10**(2.0_rk*q(2,:)-q(1,:)))/H+ gam_xxff/q10(1,:)/H)
      rhs(2,:) =  - 3.0_rk + (sv_xxaa*neqzp*neqzp*(1.0_rk/q10(2,:)- q10(2,:)/neqazp/neqazp))/H+ (gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)/q10(2,:)/H
    else
      rhs(1,:) =  - 3.0_rk +  ((-sv_xxaa*q10(1,:)+sv_aaxx* 10**(2.0_rk*q(2,:)-q(1,:)))/H+ gam_xxff/q10(1,:)/H)
      rhs(2,:) =  - 3.0_rk + (sv_xxaa*(10**(2.0_rk*q(1,:)-q(2,:))- neqzp*neqzp/neqazp/neqazp*q10(2,:)))/H+ (gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)/q10(2,:)/H
    end if

    do i = 1, params%N
      ! axions and DM in equilibrium and source term small
      if (Tp(i)>mx) then
        ! p,eq,a(T')
        peqaTp = peq(Tp(i),ma,ga)
        ! p,eq,DM(T')
        peqDMTp = peq(Tp(i),mx,gDM)
        ! dT'/dlz
        rhs(3,i) = l10*( -3.0_rk * ( rhoeqaTp + rhoeqDMTp + peqaTp + peqDMTp ) + params%gaff(i)*params%gaff(i)*drhoa/H )&
                  /(drhoeq(Tp(i),mx,gDM)+drhoeq(Tp(i),ma,ga))
      else
        rhoeqneqDM = rhoeqneq(Tp(i),mx)
        if (Tp(i)>ma) then
        ! p = peq(T')/neq(T')*n=T'*n for MB
        ! p,eq,a(T')
          peqaTp = peq(Tp(i),ma,ga)
          rhoplusp = 3.0_rk*(rhoeqaTp + peqaTp + rhoeqneqDM*q10(1,i) + Tp(i)*q10(1,i) )
          rhs(3,i) = l10*( -rhoplusp + params%gaff(i)*params%gaff(i)*drhoa/H &
                    - (rhoeqneqDM*q10(1,i)*rhs(1,i))) &
                    /(q10(1,i)*drhoeqneq( Tp(i), mx ) + drhoeq(Tp(i),ma,ga))
        else
          rhoeqneqa = rhoeqneq(Tp(i),ma)
          ! both rho/n(T')
          rhoplusp = 3.0_rk*(rhoeqneqa*q10(2,i) + rhoeqneqDM*q10(1,i) + Tp(i)*(q10(1,i)+q10(2,i)) )
          ! both rho/n(T')
          rhs(3,i) = l10*( -rhoplusp + params%gaff(i)*params%gaff(i)*drhoa/H &
                    - ( rhoeqneqa*q10(2,i)*rhs(2,i) + rhoeqneqDM*q10(1,i)*rhs(1,i))) &
                    /(q10(1,i)*drhoeqneq( Tp(i), mx ) + q10(2,i)*drhoeqneq( Tp(i), ma ))
        end if
      end if
    end do

  case ("freeze-out")
    ! DM
    ! maybe the commented out version is more stable...
    rhs(1,:) = -3.0_rk +  (sv_xxaa*neqzp*neqzp+ gam_xxff)*(1.0_rk/q10(1,:) - q10(1,:)/neqzp/neqzp)/H!((-sv_xxaa*q10(1,:)+sv_aaxx* 10**(2.0_rk*q(2,:)-q(1,:)))/H+ gam_xxff*(1.0_rk/q10(1,:) - q10(1,:)/neqzp/neqzp)/H)
    ! - 3.0_rk +  (sv_xxaa*neqzp*neqzp+ gam_xxff)*(1.0_rk/q10(1,:) - q10(1,:)/neqzp/neqzp)/H
    rhs(2,:) =  - 3.0_rk !+ (sv_aaxx*neqazp*neqazp + gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)*(1.0_rk/q10(2,:) - q10(2,:)/neqazp/neqazp)/H
    rhs(3,:) = -l10*T
  case ("freeze-in")
    rhs(1,:) =  - 3.0_rk +  (sv_aaxx*neqazp*neqazp + gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)/H/q10(1,:)
    rhs(2,:) =  - 3.0_rk + (-sv_aaxx*neqazp*neqazp/q10(2,:)+ (gam_agff + 2.0_rk*gam_afgf)*(1.0_rk/q10(2,:) - 1.0_rk/neqazp) + ffa*(neqazp/q10(2,:) - 1.0_rk))/H
    !((-sv_aaxx*neqazp*neqazp)/q10(2,:)+ (gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)*(1.0_rk/q10(2,:) - 1.0_rk/neqazp))/H
    rhs(3,:) =  - l10*T
  case ("seq-freeze-in")
    rhs(1,:) =  - 3.0_rk + (sv_aaxx*10**(2.0_rk*q(2,:)-q(1,:)) + gam_xxff/q10(1,:))/H
    rhs(2,:) =  - 3.0_rk + (-sv_aaxx*q10(2,:) + (gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)/q10(2,:))/H
    rhs(3,:) =  - l10*T
  case default
    write(*,*) "Regime ", params%regime, " not (yet) implemented."
  end select

Ynew = reshape(rhs,(/N/))
end subroutine boltzmann_logn
