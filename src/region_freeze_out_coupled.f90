subroutine region_freeze_out_cp( N, lz, Y, Ynew, params, argsint )
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
  real(kind=rk), dimension(params%N)                :: Tp, sv_aaxx,sv_xxaa, neqz,neqaz, ffa
  real(kind=rk), dimension(params%N)                :: gam_xxff, gam_agff
  real(kind=rk)                                     :: rhoeqaT,rhoeqDMT

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  ma = params%ma
  T = mx/10**lz
  ! rho,eq,a(T)
  rhoeqaT = rhoeq(T,ma,ga)
  ! rhoeq,DM(T)
  rhoeqDMT = rhoeq(T,mx,gDM)
  ! Hubble function
  H = Hub( T, rhoeqaT+rhoeqDMT )
  do i=1,params%N
    ! neq,DM(z')
    neqz(i) = neq(T, mx, gDM)
    ! neq,a(z')
    neqaz(i) = neq(T, ma, ga)
    ! HS interaction
    ! argsint%g = params%gaxx(i)
    call sigmav( T, params, argsint, "aaxx", sv_aaxx(i) )
    sv_aaxx(i) = sv_aaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    call sigmav( T, params, argsint, "xxaa", sv_xxaa(i) )
    sv_xxaa(i) = sv_xxaa(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)*params%gaxx(i)
    call gamma_r_new( T, params, argsint, "agffth", gam_agff(i) )
    gam_agff(i) = gam_agff(i)*params%gaff(i)*params%gaff(i)

    call gamma_r_new( T, params, argsint, "xxffth", gam_xxff(i) )
    !call gamma_r_new( T, params, argsint, "xxff", gam_xxff(i) )
    gam_xxff(i) = gam_xxff(i)*params%gaxx(i)*params%gaff(i)*params%gaxx(i)*params%gaff(i)

    ffa(i) = gammav(T, params, "affth")*params%gaff(i)*params%gaff(i)

  end do

  ! DM
  rhs(1,:) =  -l10*3.0_rk*q(1,:) + &
              !l10* (sv_aaxx*neqaz*neqaz*( 1.0_rk - q(1,:)*q(1,:)/neqz/neqz)/H &
              !l10* (sv_xxaa*( q(2,:)*q(2,:)*neqz*neqz/neqaz/neqaz - q(1,:)*q(1,:))/H &
              l10* ((-sv_xxaa*q(1,:)*q(1,:)+sv_aaxx*q(2,:)*q(2,:))/H &
              + gam_xxff*(1.0_rk - q(1,:)*q(1,:)/neqz/neqz)/H)
  ! axions
  rhs(2,:) =  -l10*3.0_rk*q(2,:) + &
              !l10* (sv_aaxx*neqaz*neqaz*( -1.0_rk + q(1,:)*q(1,:)/neqz/neqz)/H &
              l10* (sv_xxaa*( -q(2,:)*q(2,:)*neqz*neqz/neqaz/neqaz + q(1,:)*q(1,:))/H &
              !l10* ((sv_xxaa*q(1,:)*q(1,:)-sv_aaxx*q(2,:)*q(2,:))/H &
              + gam_agff*(1.0_rk - q(2,:)*q(2,:)/neqaz/neqaz)/H + ffa*(neqaz-q(2,:))/H)

Ynew = reshape(rhs,(/N/))
end subroutine region_freeze_out_cp
