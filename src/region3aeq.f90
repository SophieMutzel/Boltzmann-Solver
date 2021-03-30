subroutine region3a_eq( N, lz, Y, Ynew, params, argsint )
  implicit none
  real(kind=rk), intent(in)                         :: lz
  real(kind=rk), dimension(N), intent(in)           :: Y
  real(kind=rk), dimension(N), intent(inout)        :: Ynew
  type (type_params), intent(in)                    :: params
  type (type_argsint), intent(inout)                :: argsint
  integer(kind=ik), intent(in)                      :: N
  real(kind=rk), dimension(nrhs,params%N)           :: q, rhs
  real(kind=rk)                                     :: mx, T, mf, nc, s, H!, hSM, rar, gSM, gHS
  real(kind=rk)                                     :: result
  integer(kind=ik)                                  :: ier, i, neval!,nd, nr
  real(kind=rk), dimension(params%N)                :: Tprim, sv_aaxx, sv_xxaa, neqzp,neqazp
  real(kind=rk), dimension(params%N)                :: gam_agff, gam_afgf, gam_xxff

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  T = mx/10**lz
  Tprim(:) = Ta(T,params)
  ! HS interaction
  do i=1,params%N
    argsint%g = params%gaxx(i)
    call sigmav( Tprim(i), params, argsint, "aaxx", sv_aaxx(i) )
    !call sigmav( Tprim(i), params, argsint, "xxaa", sv_xxaa(i) )

    ! SM axion interaction
    !argsint%g = params%gaff(i)
    !call gamma_r_new( T, params, argsint, "agff", gam_agff(i) )
    !call gamma_r_new( T, params, argsint, "afgf", gam_afgf(i) )

    ! SM DM interaction
    ! argsint%g = params%gaxx(i)*params%gaff(i)
    !call gamma_r_new( T, params, argsint, "xxff", gam_xxff(i) )
  end do

  ! Y_x
  do i=1,params%N
    neqzp(i) = neq(Tprim(i), params%mx, gDM)
    neqazp(i) = neq(Tprim(i), params%ma, ga)
  end do
  s = ent( T, params )
  H = Hub( T, params )

  rhs(1,:) = sv_aaxx*neqazp*neqazp*(1.0_rk-q(1,:)*q(1,:)/neqzp/neqzp*s*s)/H/s
  rhs(2,:) = -sv_aaxx*neqazp*neqazp*(1.0_rk-q(1,:)*q(1,:)/neqzp/neqzp*s*s)/H/s!+(gam_agff)/H/s
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
