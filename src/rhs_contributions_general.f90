subroutine rhs_contributions_general( N, lz, Y, params, argsint, rhs )
  implicit none
  real(kind=rk), intent(in)                   :: lz
  real(kind=rk), dimension(N), intent(in)     :: Y
  real(kind=rk), dimension(:,:), intent(out)  :: rhs
  type (type_params), intent(in)              :: params
  type (type_argsint), intent(inout)          :: argsint
  integer(kind=ik), intent(in)                :: N
  real(kind=rk), dimension(nrhs,params%N)     :: q
  real(kind=rk)                               :: mx, T, mf, nc, ma, H
  real(kind=rk)                               :: result
  real(kind=rk)                               :: sv_agff, sv_afgf, sv_xxff
  integer(kind=ik)                            :: ier, i, neval,nd, nz
  real(kind=rk), dimension(params%N)          :: Tp, sv_aaxx, sv_xxaa, sv_axax,neqzp, neqazp
  real(kind=rk), dimension(params%N)          :: gam_agff, gam_afgf, gam_xxff, ffa
  real(kind=rk)                               :: rhoeqneqa,rhoeqneqDM, drhoa
  real(kind=rk)                               :: rhoeqaTp, rhoeqDMTp, peqaTp, peqDMTp,rhoplusp, rhs_a, rhs_DM

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  ma = params%ma
  T = mx/10**lz
  select case (params%regime)
  case("reannihilation")
    Tp = q(nrhs,:)
  case("freeze-out")
    Tp(:) = T
  case("freeze-in")
    Tp = T
  case("seq-freeze-in")
    Tp(:) = T
  case default
    write(*,*) "Error! Regime ", params%regime, " not (yet) implemented!"
  end select
  argsint%T = T
  ! rho,eq,a(T')
  rhoeqaTp = rhoeq(Tp(1),ma,ga)
  ! rhoeq,DM(T')
  rhoeqDMTp = rhoeq(Tp(1),mx,gDM)
  ! Hubble function
  H = Hub( T, rhoeqaTp+rhoeqDMTp )
  ! HS interaction
  do i=1,params%N
    ! DM axion interaction
    !argsint%g = params%gaxx(i)
    call sigmav( Tp(i), params, argsint, "aaxx", sv_aaxx(i) )
    call sigmav( Tp(i), params, argsint, "xxaa", sv_xxaa(i) )
    call sigmav( Tp(i), params, argsint, "axax", sv_axax(i) )
    sv_aaxx(i) = sv_aaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    sv_xxaa(i) = sv_xxaa(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    sv_axax(i) = sv_axax(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    ! SM axion interaction
    !argsint%g = params%gaff(i)
    call gamma_r_new( T, argsint, "agffth", gam_agff(i) )
    !call gamma_r_new( T, argsint, "agff", gam_agff(i) )
    gam_agff(i) = gam_agff(i)*params%gaff(i)*params%gaff(i)
    !gam_afgf(i) = 0.0_rk
    !call gamma_r_new( T, argsint, "afgfth", gam_afgf(i) )
    nz = size(params%gam_afgf,2)
    call interp_linear(nz, params%gam_afgf(1,:),params%gam_afgf(2,:),lz, gam_afgf(i))
    gam_afgf(i) = gam_afgf(i)*params%gaff(i)*params%gaff(i)
    ! SM DM interaction
    !argsint%g = params%gaxx(i)*params%gaff(i)
    call gamma_r_new( T, argsint, "xxffth", gam_xxff(i) )
    !call gamma_r_new( T, argsint, "xxff", gam_xxff(i) )
    gam_xxff(i) = gam_xxff(i)*params%gaxx(i)*params%gaff(i)*params%gaxx(i)*params%gaff(i)

    ffa(i) = gammav(T, argsint, "affth")*params%gaff(i)*params%gaff(i)
  end do
  ! Y_x
!  do i=1,params%N
!    neqzp(i) = neq(Tp(i), params%mx, gDM)
!    neqazp(i) = neq(Tp(i), params%ma, ga)
!  end do
  rhs(1,:) = lz
  rhs(2,:) = H
  rhs(3,:) = sv_xxaa*10**q(1,:)!neqzp(:)!
  rhs(4,:) = sv_aaxx*10**q(2,:)!neqazp(:)!10**q(2,:)
  rhs(5,:) = sv_axax*10**q(2,:)!neqazp(:)!10**q(2,:) ! or nx? would have to be nx for na equation
  rhs(6,:) = gam_xxff/neq(T, params%mx, gDM)!10**q(1,:)!/
  rhs(7,:) = (gam_agff + 2.0_rk*gam_afgf) /neq(T, params%ma, ga) + ffa
  !write(*,*) gam_agff, neqazp, ffa, 2.0_rk*gam_afgf
            !  +ffa*neqazp
!  rhs(1,:) = lz
!  rhs(2,:) = Tp
!  rhs(3,:) = l10*(sv_xxaa*(q(1,:)*q(1,:)))/H
!  rhs(4,:) = l10*(gam_xxff)/H
!  rhs(5,:) = l10*3.0_rk*q(1,:)
!  rhs(6,:) = l10*3.0_rk*q(2,:)
!  rhs(7,:) = l10*(gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)/H


end subroutine rhs_contributions_general
