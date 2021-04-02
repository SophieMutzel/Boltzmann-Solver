subroutine region3a_eq( N, lz, Y, Ynew, params, argsint )
  implicit none
  real(kind=rk), intent(in)                         :: lz
  real(kind=rk), dimension(N), intent(in)           :: Y
  real(kind=rk), dimension(N), intent(inout)        :: Ynew
  type (type_params), intent(in)                    :: params
  type (type_argsint), intent(inout)                :: argsint
  integer(kind=ik), intent(in)                      :: N
  real(kind=rk), dimension(nrhs,params%N)           :: q, rhs
  real(kind=rk)                                     :: mx, T, mf, nc, s, H, drhoa
  real(kind=rk)                                     :: result
  integer(kind=ik)                                  :: ier, i, neval,nd!, nr
  real(kind=rk), dimension(params%N)                :: Tprim, sv_aaxx, sv_xxaa, neqzp,neqazp
  real(kind=rk), dimension(params%N)                :: gam_agff, gam_afgf, gam_xxff
  real(kind=rk)                                     :: rhoeqDM, neqDM, neqa, rhoeqa, dTdlz, dT, ds

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  T = mx/10**lz
  s = ent( T, params )
  H = Hub( T, params )
  !Tprim(:) = Tanew(T,params,q(3,:),q(1,:)/neq(T,mx,gDM)*s*rhoeq(T,mx,gDM))
  Tprim(:) = q(3,:)
  ! HS interaction
  do i=1,params%N
    argsint%g = params%gaxx(i)
    call sigmav( Tprim(i), params, argsint, "aaxx", sv_aaxx(i) )
    !call sigmav( Tprim(i), params, argsint, "xxaa", sv_xxaa(i) )

    ! SM axion interaction
    argsint%g = params%gaff(i)
    call gamma_r_new( T, params, argsint, "agff", gam_agff(i) )
    call gamma_r_new( T, params, argsint, "afgf", gam_afgf(i) )

    ! SM DM interaction
    ! argsint%g = params%gaxx(i)*params%gaff(i)
    !call gamma_r_new( T, params, argsint, "xxff", gam_xxff(i) )
  end do

  ! Y_x
  do i=1,params%N
    neqzp(i) = neq(Tprim(i), mx, gDM)
    neqazp(i) = neq(Tprim(i), params%ma, ga)
  end do

  nd = size(params%drhoa,2)
  call interp_linear(nd, params%drhoa(1,:),params%drhoa(2,:),T, drhoa)
  rhs(1,:) = l10 * (sv_aaxx*neqazp*neqazp*(1.0_rk-q(1,:)*q(1,:)/neqzp/neqzp*s*s)/H/s)
  rhs(2,:) = l10 * (-sv_aaxx*neqazp*neqazp*(1.0_rk-q(1,:)*q(1,:)/neqzp/neqzp*s*s)/H/s + (gam_agff+2.0_rk*gam_afgf)/H/s)
!  if (Tprim(1)>mx) then
!    rhs(3,:) = -4.0_rk*q(3,:)+drhoa/H
!  else
!    rhs(3,:) = -4.0_rk*q(3,:)+mx*q(1,:)/s+drhoa/H
!  end if
  rhoeqDM = rhoeq(T,mx,gDM)
  neqDM   = neq(T,mx,gDM)
  neqa    = neq(T,params%ma,ga)
  rhoeqa  = rhoeq(T,params%ma,ga)
  dTdlz   = -T*l10
  dT      = 0.1_rk
  if ( q(3,1) > mx ) then
    rhs(3,:) =(-4.0_rk*l10*(rhoeqDM+rhoeq(T,params%ma,ga)) +l10*params%gaff*params%gaff*drhoa)/(drhoeq(T,mx,gDM)+drhoeq(T,params%ma,ga))
  else
    if (q(3,1)>params%ma) then
      ds = 2.0_rk/15.0_rk*pi*pi*T*T* heff(T, params) + s * 0.5_rk*(heff(T+dT, params)-heff(T-dT, params))/dT
      rhs(3,:) = (l10*(-4.0_rk*(rhoeq(T,params%ma,ga)+rhoeqDM/neqDM*s*q(1,:)) - mx*q(1,:)*s + drhoa*params%gaff*params%gaff) - &
                (drhoeqneq( T, mx, gDM )*dTdlz*s + dTdlz*ds*rhoeqDM/neqDM  + rhoeqDM/neqDM*s*rhs(1,:)))/(drhoeq(T,params%ma,ga))
    else
      rhs(3,:) = 0.0_rk
    end if
  end if
!if ((sv_aaxx(1)*neqazp(1) > 0.0000001_rk*H) .and. (sv_xxaa(1)*neqzp(1) > 0.0000001_rk*H) .and. ((gam_agff(1)+2.0_rk*gam_afgf(1)) < 0.0000001_rk*sv_aaxx(1)*neqazp(1)*neqazp(1))&
!    .and. gam_xxff(1) < 0.0000001_rk*sv_aaxx(1)*neqazp(1)*neqazp(1)) then
!  rhs(1,:) =  (-sv_xxaa*q(1,:)*q(1,:)+sv_aaxx* q(2,:)*q(2,:))*s/H + gam_xxff/s/H
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
