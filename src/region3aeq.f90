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
  !T= mx/lz
  argsint%T = T
!  nd = size(params%heff_HS,2)
!  ! DOFS
!  call interp_linear(nd, params%heff_HS(1,:),params%heff_HS(2,:),T, gHS)
!  call geffSM(T,params,gSM)
!  ! Axion temperature
!  nr = size(params%rhoa_rho,2)
!  call interp_linear(nr, params%rhoa_rho(1,:),params%rhoa_rho(2,:),T, rar)
!  Tprim(:) = sqrt(params%gaff*sqrt(( gSM*gSM/gHS *rar )))*T
  Tprim(:) = sqrt(params%gaff) * Ta(T,params)
  ! HS interaction
  do i=1,params%N
    argsint%g = params%gaxx(i)
    call sigmav( Tprim(i), params, argsint, "aaxx", sv_aaxx(i) )
    call sigmav( Tprim(i), params, argsint, "xxaa", sv_xxaa(i) )
    !call sigmav( T, params, argsint, "aaxx", sv_aaxx(i) )
    !call sigmav( T, params, argsint, "xxaa", sv_xxaa(i) )
    ! SM axion interaction
    argsint%g = params%gaff(i)
    call gamma_r_new( T, params, argsint, "agff", gam_agff(i) )
    call gamma_r_new( T, params, argsint, "afgf", gam_afgf(i) )

    ! SM DM interaction
     argsint%g = params%gaxx(i)*params%gaff(i)
    call gamma_r_new( T, params, argsint, "xxff", gam_xxff(i) )
  end do

  ! Y_x
  do i=1,params%N
    neqzp(i) = neq(sqrt(params%gaff(i))*Ta(T,params), params%mx, gDM)
    neqazp(i) = neq(sqrt(params%gaff(i))*Ta(T,params), params%ma, ga)
  end do
  s = ent( T, params )
  H = Hub( T, params )
  !rhs(1,:) = params%gaxx**4 * (sv_xxaa*neqzp*neqzp/s&
  !          - sv_xxaa*q(1,:)*q(1,:)*s)/H!+&
            !params%gaff*params%gaff*params%gaxx*params%gaxx*&
            !gam_xxff/s/H
  !rhs(1,:) = params%gaxx**4 * (-sv_aaxx*neqazp*neqazp/s&
  !          + sv_aaxx*q(2,:)*q(2,:)*s)/H!+&
            !params%gaff*params%gaff*params%gaxx*params%gaxx*&
            !gam_xxff/s/H
  ! Y_a
  !rhs(2,:) = 1.0_rk/H * params%gaxx**4 * (-sv_xxaa*neqzp*neqzp/s &
  !          + sv_xxaa*q(1,:)*q(1,:)*s )&
  !          + params%gaff*params%gaff*(gam_agff + gam_afgf)/s/H
  !rhs(2,:) = params%gaxx**4 * (-sv_aaxx*q(2,:)*q(2,:)*s &
  !          + sv_aaxx*neqazp*neqazp/s )/H&
  !          + params%gaff*params%gaff*(gam_agff + gam_afgf)/s/H


  !if ((sv_aaxx(1)*neqazp(1) > 0.01_rk*H) .and. (sv_xxaa(1)*neqzp(1) > 0.01_rk*H) .and. ((gam_agff(1)+gam_afgf(1)) < 0.01_rk*sv_aaxx(1)*neqazp(1)*neqazp(1))&
  !    .and. gam_xxff(1) < 0.0001_rk*sv_aaxx(1)*neqazp(1)*neqazp(1)) then
    !write(*,*) lz
!  if (lz <-0.3_rk) then
!    rhs(1,:) =  sv_aaxx*neqazp*neqazp*(1.0_rk&
!              - q(1,:)*q(1,:)/neqzp/neqzp*s*s )/H/s+ gam_xxff/H/s
!    !rhs(2,:) =  -sv_aaxx*(1.0_rk&
!    !          - q(1,:)*q(1,:)/neqzp/neqzp*s*s)/H/s
!    rhs(2,:) =  sv_xxaa*neqzp*neqzp*(1.0_rk&
!              - q(2,:)*q(2,:)/neqazp/neqazp*s*s)/H/s + (gam_agff + gam_afgf)/H/s
!  else
!    rhs(1,:) =  sv_xxaa*(-q(1,:)*q(1,:)&
!              !+ sv_aaxx*neq(T,params%ma,ga)*neq(T,params%ma,ga)/s)/H
!            + q(2,:)*q(2,:)*neqzp*neqzp/neqazp/neqazp)*s/H + gam_xxff/s/H
!    rhs(2,:) =  sv_xxaa*(q(1,:)*q(1,:)&
!              !- sv_aaxx*neq(T,params%ma,ga)*neq(T,params%ma,ga)/s)/H
!            - q(2,:)*q(2,:)*neqzp*neqzp/neqazp/neqazp)*s/H &
!            + (gam_agff + gam_afgf)/s/H
!  end if
if ((sv_aaxx(1)*neqazp(1) > 0.01_rk*H) .and. (sv_xxaa(1)*neqzp(1) > 0.01_rk*H) .and. ((gam_agff(1)+gam_afgf(1)) < 0.01_rk*sv_aaxx(1)*neqazp(1)*neqazp(1))&
    .and. gam_xxff(1) < 0.0001_rk*sv_aaxx(1)*neqazp(1)*neqazp(1)) then
  rhs(1,:) =  (-sv_xxaa*q(1,:)*q(1,:)+sv_aaxx* q(2,:)*q(2,:))*s/H + gam_xxff/s/H
  rhs(2,:) =  sv_xxaa*neqzp*neqzp*(1.0_rk - q(2,:)*q(2,:)/neqazp/neqazp*s*s)/H/s+ (gam_agff + gam_afgf)/H/s
else
  rhs(1,:) =  sv_xxaa*(-q(1,:)*q(1,:)&
            !+ sv_aaxx*neq(T,params%ma,ga)*neq(T,params%ma,ga)/s)/H
            + q(2,:)*q(2,:)*neqzp*neqzp/neqazp/neqazp)*s/H + gam_xxff/s/H
  rhs(2,:) =  sv_xxaa*(q(1,:)*q(1,:)&
            !- sv_aaxx*neq(T,params%ma,ga)*neq(T,params%ma,ga)/s)/H
            - q(2,:)*q(2,:)*neqzp*neqzp/neqazp/neqazp)*s/H &
            + (gam_agff + gam_afgf)/s/H
end if
Ynew = reshape(rhs,(/N/))
end subroutine region3a_eq
