subroutine region3a_eq( N, lz, Y, Ynew, params, argsint )
  implicit none
  real(kind=rk), intent(in)                         :: lz
  real(kind=rk), dimension(N), intent(in)           :: Y
  real(kind=rk), dimension(N), intent(inout)        :: Ynew
  type (type_params), intent(in)                    :: params
  type (type_argsint), intent(inout)                :: argsint
  integer(kind=ik), intent(in)                      :: N
  real(kind=rk), dimension(nrhs,params%N)           :: q, rhs
  real(kind=rk)                                     :: mx, ma, T, mf, nc, s, H, drhoa
  real(kind=rk)                                     :: result
  integer(kind=ik)                                  :: ier, i, neval,nd!, nr
  real(kind=rk), dimension(params%N)                :: Tp, sv_aaxx, sv_xxaa, neqzp,neqazp
  real(kind=rk), dimension(params%N)                :: gam_agff, gam_afgf, gam_xxff
  real(kind=rk)                                     :: rhoeqDMT, neqDM, neqa, rhoeqaT, dT, ds
  real(kind=rk)                                     :: rhoeqaTp, rhoeqDMTp, peqaTp, peqDMTp,rhoplusp

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  ma = params%ma
  T = mx/10**lz
  s = ent( T, params )
  Tp(:) = q(3,:)!Ta(T,params)
  ! rho,eq,a(T')
  rhoeqaTp = rhoeq(Tp(1),ma,ga)
  ! rhoeq,DM(T')
  rhoeqDMTp = rhoeq(Tp(1),mx,gDM)
  ! Hubble function
  H = Hub( T, rhoeqaTp+rhoeqDMTp )
  !Tprim(:) = Tanew(T,params,q(3,:),q(1,:)/neq(T,mx,gDM)*s*rhoeq(T,mx,gDM))
  do i=1,params%N
    ! HS interaction
    argsint%g = params%gaxx(i)
    call sigmav( Tp(i), params, argsint, "aaxx", sv_aaxx(i) )
    call sigmav( Tp(i), params, argsint, "xxaa", sv_xxaa(i) )

    ! SM axion interaction
    argsint%g = params%gaff(i)
    call gamma_r_new( T, params, argsint, "agff", gam_agff(i) )
    call gamma_r_new( T, params, argsint, "afgf", gam_afgf(i) )

    ! SM DM interaction
    ! argsint%g = params%gaxx(i)*params%gaff(i)
    !call gamma_r_new( T, params, argsint, "xxff", gam_xxff(i) )
  end do
  do i=1,params%N
    ! neq,DM(z')
    neqzp(i) = neq(Tp(i), mx, gDM)
    ! neq,a(z')
    neqazp(i) = neq(Tp(i), params%ma, ga)
  end do
  ! derivative of entropy with respect to T
  ds = 2.0_rk/15.0_rk*pi*pi*T*T*heff(T,params)+2.0_rk/45.0_rk*pi*pi*T*T*T*0.5_rk*(heff(T+dT, params)-heff(T-dT, params))/dT
  ! DM
  !rhs(1,:) =  (sv_aaxx*neqazp*neqazp*(1.0_rk-q(1,:)*q(1,:)/neqzp/neqzp*s*s)/H/s)
  ! axion
  !rhs(2,:) =  sv_aaxx*neqazp*neqazp*(q(1,:)*q(1,:)/neqzp/neqzp*s*s-1.0_rk)/H/s + (gam_agff + 2.0_rk*gam_afgf)/H/s
  !write(*,*)  q(1,:)*q(1,:)/neqzp/neqzp*s*s, q(2,:)*q(2,:)/neqazp/neqazp*s*s
  !rhs(2,:) =  sv_aaxx*neqazp*neqazp*(1.0_rk - q(2,:)*q(2,:)/neqazp/neqazp*s*s)/H/s+ (gam_agff + 2.0_rk*gam_afgf)/H/s
  rhs(1,:) =  l10* (-sv_xxaa*q(1,:)*q(1,:)+sv_aaxx* q(2,:)*q(2,:))*s/H !+ gam_xxff/s/H
  rhs(2,:) =  l10*(sv_xxaa*neqzp*neqzp*(1.0_rk - q(2,:)*q(2,:)/neqazp/neqazp*s*s)/H/s+ (gam_agff + 2.0_rk*gam_afgf)/H/s)
  ! collision term for energy transfer
  nd = size(params%drhoa,2)
  call interp_linear(nd, params%drhoa(1,:),params%drhoa(2,:),T, drhoa)
  ! equilibrium values
  ! rho,eq,DM(T)
  rhoeqDMT = rhoeq(T,mx,gDM)
  ! rho,eq,a(T)
  rhoeqaT  = rhoeq(T,ma,ga)
  ! neq,DM(T)
  neqDM    = neq(T,mx,gDM)
  ! neq,a(T)
  neqa     = neq(T,ma,ga)
  ! delta T for finite difference in heff
  dT      = 1e-8_rk
  do i = 1, params%N
    ! rho,eq,a(T')
    rhoeqaTp = rhoeq(Tp(i),ma,ga)
    ! rhoeq,DM(T')
    rhoeqDMTp = rhoeq(Tp(i),mx,gDM)
    ! p,eq,a(T')
    peqaTp = peq(Tp(i),ma,ga)
    ! p,eq,DM(T')
    peqDMTp = peq(Tp(i),mx,gDM)
    ! axions and DM in equilibrium and source term small
    !if ((sv_aaxx(i)*neqazp(i) > 0.001_rk*H) .and. (sv_xxaa(i)*neqzp(i) > 0.001_rk*H) &
    !    .and. ((gam_agff(i)+2.0_rk*gam_afgf(i)) < 0.001_rk*sv_aaxx(i)*neqazp(i)*neqazp(i))) then
    if (lz<0.0_rk) then
      ! dT'/dlz
      rhs(3,i) = l10*( -3.0_rk * ( rhoeqaTp + rhoeqDMTp + peqaTp + peqDMTp ) - params%gaff(i)*params%gaff(i)*drhoa/H )&
                /(drhoeq(Tp(i),mx,gDM)+drhoeq(Tp(i),ma,ga))
    else
      ! rho = rhoeq(T')/neq(T')*s*Y
      ! p = peq(T')/neq(T')*s*Y=T'*s*Y for MB
      ! both rho/n(T')
      rhoplusp = 3.0_rk*s*(rhoeqaTp/neqazp(i)*q(2,i) + rhoeqDMTp/neqzp(i)*q(1,i) + Tp(i)*(q(1,i)+q(2,i)) )
      ! DM: rho/n(T'), a: rho/n(T)
      !rhoplusp = 3.0_rk*s*(rhoeqaT/neqa*q(2,i) + rhoeqDMTp/neqzp(i)*q(1,i) + Tp(i)*q(1,i)+T*q(2,i) )
      ! DM: rho/n(T), a: rho/n(T')
      !rhoplusp = 3.0_rk*s*(rhoeqaTp/neqazp(i)*q(2,i) + rhoeqDMT/neqDM*q(1,i) + T*q(1,i)+Tp(i)*q(2,i) )
      ! a: equilibrium distribution, DM: rho/n(T')
      !rhoplusp = 3.0_rk*(rhoeqaTp + s*rhoeqDMTp/neqzp(i)*q(1,i) + Tp(i)*(q(1,i)*s+neqazp(i)) )
      ! a: equilibrium distribution, DM: rho/n(T)
      !rhoplusp = 3.0_rk*(rhoeqaTp +peqaTp+ s*rhoeqDMT/neqDM*q(1,i) + T*q(1,i)*s )
      ! dT'/dlz
      ! both rho/n(T')
      rhs(3,i) = (l10*( -rhoplusp - params%gaff(i)*params%gaff(i)*drhoa/H) &
                  -s*(rhoeqaTp/neqazp(i)*rhs(2,i) + rhoeqDMTp/neqzp(i)*rhs(1,i)) &
                  +T*l10*ds*(rhoeqaTp/neqazp(i)*q(2,i) + rhoeqDMTp/neqzp(i)*q(1,i)))&
                  /(s*(q(1,i)*drhoeqneq( Tp(i), mx) + q(2,i)*drhoeqneq( Tp(i), ma )))
      ! a: equilibrium distribution, DM: rho/n(T')
!      rhs(3,i) = (l10*( -rhoplusp - params%gaff(i)*params%gaff(i)*drhoa/H) &
!            -s*(rhoeqDMTp/neqzp(i)*rhs(1,i)) &
!            +T*l10*ds*(rhoeqDMTp/neqzp(i)*q(1,i)))&
!            /(s*(q(1,i)*drhoeqneq( Tp(i), mx, gDM )) +drhoeq(Tp(i),ma,ga))
      ! DM: rho/n(T'), a: rho/n(T)
!      rhs(3,i) = (l10*( -rhoplusp - params%gaff(i)*params%gaff(i)*drhoa/H) &
!            -s*(rhoeqDMTp/neqzp(i)*rhs(1,i)) &
!            +T*l10*ds*(rhoeqDMTp/neqzp(i)*q(1,i))+&
!            T*l10*drhoeqneq( T, ma, ga )*s*q(2,i)+T*l10*rhoeqaT/neqa*q(2,i)*ds&
!                          -s*rhoeqaT/neqa*rhs(2,i))&
!            /(s*(q(1,i)*drhoeqneq( Tp(i), mx, gDM )))

      ! DM: rho/n(T), a: rho/n(T')
!      rhs(3,i) = (l10*( -rhoplusp - params%gaff(i)*params%gaff(i)*drhoa/H) &
!            -s*(rhoeqaTp/neqazp(i)*rhs(2,i)) &
!            +T*l10*ds*(rhoeqaTp/neqazp(i)*q(2,i))+&
!            T*l10*drhoeqneq( T, mx, gDM )*s*q(1,i)+T*l10*rhoeqDMT/neqDM*q(1,i)*ds&
!                          -s*rhoeqDMT/neqDM*rhs(1,i))&
!            /(s*(q(2,i)*drhoeqneq( Tp(i), ma, ga )))
      ! DM: rho/n(T), a: equilibrium
!      rhs(3,i) = (T*l10*drhoeqneq( T, mx, gDM )*s*q(1,i)+T*l10*rhoeqDMT/neqDM*q(1,i)*ds&
!              -rhoeqDMT/neqDM*rhs(1,i)*s-l10*rhoplusp-l10*params%gaff(i)*params%gaff(i)*drhoa/H)&
!              /drhoeq(Tp(i),ma,ga)
    end if
  end do
!  do i = 1, params%N
!    ! rho,eq,a(T') and p,eq,a(T')
!    rhoeqaTp = rhoeq(Tp(i),ma,ga)
!    ! rho,eq,DM(T') and p,eq,DM(T')
!    rhoeqDMTp = rhoeq(Tp(i),mx,gDM)
!    ! dT'/dlz
!    if ( q(3,i) > mx ) then
!      peqaTp = peq(Tp(i),ma,ga)
!      peqDMTp = peq(Tp(i),mx,gDM)
!      rhs(3,i) = l10*( -3.0_rk * ( rhoeqaTp + rhoeqDMTp + peqaTp + peqDMTp ) - params%gaff(i)*params%gaff(i)*drhoa/H )&
!                /(drhoeq(Tp(i),mx,gDM)+drhoeq(Tp(i),ma,ga))
!    else
!      if ( q(3,i) > params%ma ) then
!        ds = 2.0_rk/15.0_rk*pi*pi*T*T*heff(T,params)+2.0_rk/45.0_rk*pi*pi*T*T*T*0.5_rk*(heff(T+dT, params)-heff(T-dT, params))/dT
!        rhoplusp = 4.0_rk*rhoeqaTp + 3.0_rk*rhoeqDMT/neqDM*s*q(1,i)
!        rhs(3,i) = (T*l10*drhoeqneq( T, mx, gDM )*s*q(1,i)+T*l10*rhoeqDMT/neqDM*q(1,i)*ds&
!                  -rhoeqDMT/neqDM*rhs(1,i)-l10*rhoplusp+l10*params%gaff(i)*params%gaff(i)*drhoa/H)&
!                  /drhoeq(Tp(i),ma,ga)
!
!      else
!        if ( q(3,i) > params%ma/5.0_rk ) then
!          ds = 2.0_rk/15.0_rk*pi*pi*T*T*heff(T,params)+2.0_rk/45.0_rk*pi*pi*T*T*T*0.5_rk*(heff(T+dT, params)-heff(T-dT, params))/dT
!          rhoplusp = 3.0_rk*(rhoeqaTp+rhoeqDMT/neqDM*s*q(1,i))
!          rhs(3,i) = (T*l10*drhoeqneq( T, mx, gDM )*s*q(1,i)+T*l10*rhoeqDMT/neqDM*q(1,i)*ds&
!                  -rhoeqDMT/neqDM*rhs(1,i)-l10*rhoplusp+l10*params%gaff(i)*params%gaff(i)*drhoa/H)&
!                  /drhoeq(Tp(i),ma,ga)
!        else
!          rhs(3,i) = 0.0_rk
!        end if
!      end if
!    end if
!  end do
!  if ( q(3,1) > mx ) then
    !rhs(3,:) =(-4.0_rk*l10*(rhoeqDM+rhoeq(T,params%ma,ga)) +l10*params%gaff*params%gaff*drhoa)/(drhoeq(Tprim(1),mx,gDM)+drhoeq(Tprim(1),params%ma,ga))
!  else
!    if ( q(3,1) > params%ma ) then
!      ds = 2.0_rk/15.0_rk*pi*pi*T*T* heff(T, params) + s * 0.5_rk*(heff(T+dT, params)-heff(T-dT, params))/dT
!      rhs(3,:) = (l10*(-4.0_rk*(rhoeq(T,params%ma,ga)+rhoeqDM/neqDM*s*q(1,:)) - mx*q(1,:)*s + drhoa*params%gaff*params%gaff) - &
!                (drhoeqneq( T, mx, gDM )*dTdlz*s + dTdlz*ds*rhoeqDM/neqDM  + rhoeqDM/neqDM*s*rhs(1,:)))/(drhoeq(T,params%ma,ga))
!    else
!      rhs(3,:) = 0.0_rk
!    end if
!  end if
!if ((sv_aaxx(1)*neqazp(1) > 0.0000001_rk*H) .and. (sv_xxaa(1)*neqzp(1) > 0.0000001_rk*H) .and. ((gam_agff(1)+2.0_rk*gam_afgf(1)) < 0.0000001_rk*sv_aaxx(1)*neqazp(1)*neqazp(1))&
!    .and. gam_xxff(1) < 0.0000001_rk*sv_aaxx(1)*neqazp(1)*neqazp(1)) then
!  rhs(1,:) =  (-sv_xxaa*q(1,:)*q(1,:)+sv_aaxx* q(2,:)*q(2,:))*s/H !+ gam_xxff/s/H
!  rhs(2,:) =  sv_xxaa*neqzp*neqzp*(1.0_rk - q(2,:)*q(2,:)/neqazp/neqazp*s*s)/H/s+ (gam_agff + 2.0_rk*gam_afgf)/H/s
!else
!  rhs(1,:) =  sv_xxaa*(-q(1,:)*q(1,:)&
!            !+ sv_aaxx*neq(T,params%ma,ga)*neq(T,params%ma,ga)/s)/H
!            + q(2,:)*q(2,:)*neqzp*neqzp/neqazp/neqazp)*s/H + gam_xxff/s/H
!  rhs(2,:) =  sv_xxaa*(q(1,:)*q(1,:)&
!            !- sv_aaxx*neq(T,params%ma,ga)*neq(T,params%ma,ga)/s)/H
!            - q(2,:)*q(2,:)*neqzp*neqzp/neqazp/neqazp)*s/H &
!            + (gam_agff + 2.0_rk*gam_afgf)/s/H
!end if
Ynew = reshape(rhs,(/N/))
end subroutine region3a_eq
