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
  integer(kind=ik)                                  :: ier, i, neval,nd
  real(kind=rk), dimension(params%N)                :: Tp, sv_aaxx, sv_xxaa, neqzp,neqazp
  real(kind=rk), dimension(params%N)                :: gam_agff, gam_afgf, gam_xxff, ffa
  real(kind=rk)                                     :: rhoeqneqa,rhoeqneqDM
  real(kind=rk)                                     :: rhoeqaTp, rhoeqDMTp, peqaTp, peqDMTp,rhoplusp

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  ma = params%ma
  T = mx/10**lz
  Tp(:) = q(3,:)
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
    ! argsint%g = params%gaxx(i)
    call sigmav( Tp(i), params, argsint, "aaxx", sv_aaxx(i) )
    call sigmav( Tp(i), params, argsint, "xxaa", sv_xxaa(i) )
    sv_aaxx(i) = sv_aaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    sv_xxaa(i) = sv_xxaa(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    ! SM axion interaction
    ! argsint%g = params%gaff(i)
    call gamma_r_new( T, params, argsint, "agffth", gam_agff(i) )
    !call gamma_r_new( T, params, argsint, "agff", gam_agff(i) )
    gam_agff(i) = gam_agff(i)*params%gaff(i)*params%gaff(i)
    gam_afgf(i) = 0.0_rk
    !call gamma_r_new( T, params, argsint, "afgf", gam_afgf(i) )
    ! SM DM interaction
    ! argsint%g = params%gaxx(i)*params%gaff(i)
    call gamma_r_new( T, params, argsint, "xxffth", gam_xxff(i) )
    !call gamma_r_new( T, params, argsint, "xxff", gam_xxff(i) )
    gam_xxff(i) = gam_xxff(i)*params%gaxx(i)*params%gaff(i)*params%gaxx(i)*params%gaff(i)
    ! inverse decay ff->a
    ffa(i) = gammav(T, params, "affth")*params%gaff(i)*params%gaff(i)
    !ffa(i)=0.0_rk
  end do

!if (lz<-0.2_rk) then
!if (argsint%helper .eqv. .false.) then
!if (ma<Tp(1)) then
!!!!  ! DM
!  rhs(1,:) =  -l10*3.0_rk*q(1,:) + l10* ((-sv_xxaa*q(1,:)*q(1,:)+sv_aaxx* q(2,:)*q(2,:))/H+ gam_xxff/H)
!  !rhs(1,:) =  -l10*3.0_rk*q(1,:) + l10* (sv_aaxx*neqazp*neqazp*( 1.0_rk - q(1,:)*q(1,:)/neqzp/neqzp)/H+ gam_xxff/H)
!  !rhs(2,:) =  -l10*3.0_rk*q(2,:) + l10*(sv_aaxx*neqazp*neqazp*(-1.0_rk+q(1,:)*q(1,:)/neqzp/neqzp)/H + (gam_agff + 2.0_rk*gam_afgf+ffa*neqazp)/H)
!  rhs(2,:) =  -l10*3.0_rk*q(2,:) + l10*((sv_xxaa*neqzp*neqzp - sv_aaxx*q(2,:)*q(2,:))/H+ (gam_agff + 2.0_rk*gam_afgf + ffa)/H)
!else
  rhs(1,:) =  -l10*3.0_rk*q(1,:) + l10* ((-sv_xxaa*q(1,:)*q(1,:)+sv_aaxx* q(2,:)*q(2,:))/H+ gam_xxff/H)
  rhs(2,:) =  -l10*3.0_rk*q(2,:) + l10*(sv_xxaa*(q(1,:)*q(1,:)-neqzp*neqzp/neqazp/neqazp* q(2,:)*q(2,:))/H+ (gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)/H)
!end if
  !rhs(1,:) =  -l10*3.0_rk*q(1,:) + l10* (sv_xxaa*(-q(1,:)*q(1,:)+neqzp*neqzp/neqazp/neqazp* q(2,:)*q(2,:))/H)!+ gam_xxff/H)
  !rhs(1,:) =  -l10*3.0_rk*q(1,:) + l10* ((sv_aaxx*neqazp*neqazp-sv_aaxx* q(2,:)*q(2,:))/H)! + gam_xxff/H)
  !rhs(1,:) =  -l10*3.0_rk*q(1,:) + l10* (sv_aaxx*neqazp*neqazp*( 1.0_rk - q(1,:)*q(1,:)/neqzp/neqzp)/H)
  !rhs(1,:) =  -l10*3.0_rk*q(1,:) + l10* (sv_xxaa*neqzp*neqzp*( -1.0_rk + q(2,:)*q(2,:)/neqazp/neqazp)/H)
  ! axions
  !if ((sv_aaxx(1)*neqazp(1) > 0.001_rk*H) .and. (sv_xxaa(1)*neqzp(1) > 0.001_rk*H) &
  !    .and. ((gam_agff(1)+2.0_rk*gam_afgf(1)+ffa(1)) < 0.001_rk*sv_aaxx(1)*neqazp(1)*neqazp(1))) then
  !if (lz<-0.5_rk) then
    !rhs(2,:) =  -l10*3.0_rk*q(2,:) + l10*((sv_xxaa*neqzp*neqzp - sv_aaxx*q(2,:)*q(2,:))/H)!)+ (gam_agff + 2.0_rk*gam_afgf + ffa)/H)
  !rhs(2,:) =  -l10*3.0_rk*q(2,:) + l10*(sv_aaxx*neqazp*neqazp*(-1.0_rk+q(1,:)*q(1,:)/neqzp/neqzp)/H )!+ (gam_agff + 2.0_rk*gam_afgf)/H)
  !rhs(2,:) =  -l10*3.0_rk*q(2,:) + l10* (sv_xxaa*neqzp*neqzp*( 1.0_rk - q(2,:)*q(2,:)/neqazp/neqazp)/H)
  !else
    !rhs(2,:) =  -l10*3.0_rk*q(2,:) + l10*((sv_xxaa*q(1,:)*q(1,:) - sv_aaxx*q(2,:)*q(2,:))/H)!+ (gam_agff + 2.0_rk*gam_afgf + ffa)/H)
  !end if
  ! collision term for energy transfer
  nd = size(params%drhoa,2)
  call interp_linear(nd, params%drhoa(1,:),params%drhoa(2,:),log10(T), drhoa)
  ! inverse decays ff->a
  drhoa = drhoa + drho_decay(T, params%ma, "ffath")

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
        rhoplusp = 3.0_rk*(rhoeqaTp + peqaTp + rhoeqneqDM*q(1,i) + Tp(i)*q(1,i) )
        rhs(3,i) = (l10*( -rhoplusp + params%gaff(i)*params%gaff(i)*drhoa/H) &
                  - (rhoeqneqDM*rhs(1,i))) &
                  /(q(1,i)*drhoeqneq( Tp(i), mx ) + drhoeq(Tp(i),ma,ga))
      else
        rhoeqneqa = rhoeqneq(Tp(i),ma)
        ! both rho/n(T')
        rhoplusp = 3.0_rk*(rhoeqneqa*q(2,i) + rhoeqneqDM*q(1,i) + Tp(i)*(q(1,i)+q(2,i)) )
        ! both rho/n(T')
        rhs(3,i) = (l10*( -rhoplusp + params%gaff(i)*params%gaff(i)*drhoa/H) &
                  - ( rhoeqneqa*rhs(2,i) + rhoeqneqDM*rhs(1,i))) &
                  /(q(1,i)*drhoeqneq( Tp(i), mx ) + q(2,i)*drhoeqneq( Tp(i), ma ))
      !rhoplusp = 3.0_rk*(rhoeqneq(T,ma)*q(2,i) + rhoeqneqDM*q(1,i) + Tp(i)*q(1,i)+T*q(2,i) )

!        rhs(3,i) = (l10*( -rhoplusp - params%gaff(i)*params%gaff(i)*drhoa/H) &
!            -rhoeqneqDM*rhs(1,i)+T*l10*drhoeqneq( T, ma)*q(2,i)&
!                          -rhoeqneq(T,ma)*rhs(2,i))&
!            /(q(1,i)*drhoeqneq( Tp(i), mx))
      end if
    end if
  end do
Ynew = reshape(rhs,(/N/))
end subroutine region3a_in_n
