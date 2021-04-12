subroutine region3a_in_n( N, lz, Y, Ynew, params, argsint )
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
  real(kind=rk)                                     :: rhoeqneqa,rhoeqneqDM
  real(kind=rk)                                     :: rhoeqaTp, rhoeqDMTp, peqaTp, peqDMTp,rhoplusp

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  ma = params%ma
  T = mx/10**lz
  H = Hub( T, params )
  Tp(:) = q(3,:)
  do i=1,params%N
    ! neq,DM(z')
    neqzp(1) = neq(Tp(1), mx, gDM)
    ! neq,a(z')
    neqazp(1) = neq(Tp(1), params%ma, ga)
    ! HS interaction
    argsint%g = params%gaxx(i)
    call sigmav( Tp(i), params, argsint, "aaxx", sv_aaxx(i) )
    call sigmav( Tp(i), params, argsint, "xxaa", sv_xxaa(i) )
    !if(lz>1.5_rk) write(*,*) "1"
    !sv_xxaa(i) = sv_aaxx(i)*neqazp(i)*neqazp(i)/neqzp(i)/neqzp(i)
    ! SM axion interaction
    argsint%g = params%gaff(i)
    call gamma_r_new( T, params, argsint, "agff", gam_agff(i) )
    call gamma_r_new( T, params, argsint, "afgf", gam_afgf(i) )
    ! SM DM interaction
    ! argsint%g = params%gaxx(i)*params%gaff(i)
    !call gamma_r_new( T, params, argsint, "xxff", gam_xxff(i) )
  end do
  ! DM
  rhs(1,:) =  l10*3.0_rk*q(1,:) + l10* (-sv_xxaa*q(1,:)*q(1,:)+sv_aaxx* q(2,:)*q(2,:))/H !+ gam_xxff/s/H
  ! axions
  if ( lz<1.5_rk ) then
    rhs(2,:) =  l10*3.0_rk*q(2,:) + l10*(sv_xxaa*neqzp*neqzp*(1.0_rk - q(2,:)*q(2,:)/neqazp/neqazp)/H+ (gam_agff + 2.0_rk*gam_afgf)/H)
  else
    rhs(2,:) =  l10*3.0_rk*q(2,:) + l10*((sv_xxaa*q(1,:)*q(1,:) - sv_aaxx*q(2,:)*q(2,:))/H + (gam_agff + 2.0_rk*gam_afgf)/H)
  end if
  ! collision term for energy transfer
  nd = size(params%drhoa,2)
  call interp_linear(nd, params%drhoa(1,:),params%drhoa(2,:),T, drhoa)
  ! equilibrium values
  ! rho,eq,DM(T)
  !rhoeqDMT = rhoeq(T,mx,gDM)
  ! rho,eq,a(T)
  !rhoeqaT  = rhoeq(T,ma,ga)
  ! neq,DM(T)
  !neqDM    = neq(T,mx,gDM)
  ! neq,a(T)
  !neqa     = neq(T,ma,ga)
  ! delta T for finite difference in heff
  !dT      = 1e-8_rk
  do i = 1, params%N
    ! axions and DM in equilibrium and source term small
    !if ((sv_aaxx(i)*neqazp(i) > 0.001_rk*H) .and. (sv_xxaa(i)*neqzp(i) > 0.001_rk*H) &
    !    .and. ((gam_agff(i)+2.0_rk*gam_afgf(i)) < 0.001_rk*sv_aaxx(i)*neqazp(i)*neqazp(i))) then
    if (lz<0.0_rk) then
      ! rho,eq,a(T')
      rhoeqaTp = rhoeq(Tp(i),ma,ga)
      ! rhoeq,DM(T')
      rhoeqDMTp = rhoeq(Tp(i),mx,gDM)
      ! p,eq,a(T')
      peqaTp = peq(Tp(i),ma,ga)
      ! p,eq,DM(T')
      peqDMTp = peq(Tp(i),mx,gDM)
      ! dT'/dlz
      rhs(3,i) = l10*( -3.0_rk * ( rhoeqaTp + rhoeqDMTp + peqaTp + peqDMTp ) - params%gaff(i)*params%gaff(i)*drhoa/H )&
                /(drhoeq(Tp(i),mx,gDM)+drhoeq(Tp(i),ma,ga))
    else
      ! derivative of entropy with respect to T
      !ds = 2.0_rk/15.0_rk*pi*pi*T*T*heff(T,params)+2.0_rk/45.0_rk*pi*pi*T*T*T*0.5_rk*(heff(T+dT, params)-heff(T-dT, params))/dT
      ! rho = rhoeq(T')/neq(T')*s*Y
      ! p = peq(T')/neq(T')*s*Y=T'*s*Y for MB
      rhoeqneqa = rhoeqneq(Tp(i),ma)
      rhoeqneqDM = rhoeqneq(Tp(i),mx)
      ! both rho/n(T')
      rhoplusp = 3.0_rk*(rhoeqneqa*q(2,i) + rhoeqneqDM*q(1,i) + Tp(i)*(q(1,i)+q(2,i)) )
      ! both rho/n(T')
      rhs(3,i) = (l10*( -rhoplusp - params%gaff(i)*params%gaff(i)*drhoa/H) &
                  - ( rhoeqneqa*rhs(2,i) + rhoeqneqDM*rhs(1,i))) &
                  /(q(1,i)*drhoeqneq( Tp(i), mx ) + q(2,i)*drhoeqneq( Tp(i), ma ))
      ! a: equilibrium distribution, DM: rho/n(T')
!      rhs(3,i) = (l10*( -rhoplusp - params%gaff(i)*params%gaff(i)*drhoa/H) &
!            -s*(rhoeqDMTp/neqzp(i)*rhs(1,i)) &
!            +T*l10*ds*(rhoeqDMTp/neqzp(i)*q(1,i)))&
!            /(s*(q(1,i)*drhoeqneq( Tp(i), mx, gDM )) +drhoeq(Tp(i),ma,ga))
    end if
  end do
Ynew = reshape(rhs,(/N/))
end subroutine region3a_in_n
