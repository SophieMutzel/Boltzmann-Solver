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
  real(kind=rk)                                     :: result, gam_agff, gam_afgf, gam_xxff
  integer(kind=ik)                                  :: ier, i, neval!,nd, nr
  real(kind=rk), dimension(params%N)                :: Tprim, sv_aaxx, sv_xxaa, neqzp,neqazp

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
    call sigmav( Tprim(i), params, argsint, "aaxx", sv_aaxx(i) )
    call sigmav( Tprim(i), params, argsint, "xxaa", sv_xxaa(i) )
  end do
  ! SM axion interaction
  call gamma_r_new( T, params, argsint, "agff", gam_agff )
  call gamma_r_new( T, params, argsint, "afgf", gam_afgf )

  ! SM DM interaction
  call gamma_r_new( T, params, argsint, "xxff", gam_xxff )
  ! Y_x
  do i=1,params%N
    neqzp(i) = neq(sqrt(params%gaff(i))*Ta(T,params), params%mx, gDM)
    neqazp(i) = neq(sqrt(params%gaff(i))*Ta(T,params), params%ma, ga)
  end do
  s = ent( T, params )
  H = Hub( T, params )
  rhs(1,:) = 1.0_rk/H *params%gaxx**4 * (sv_xxaa*neqzp*neqzp/s&
            - sv_xxaa*q(1,:)*q(1,:)*s)!+&
            !params%gaff*params%gaff*params%gaxx*params%gaxx*&
            !gam_xxff/s/H
  ! Y_a
  rhs(2,:) = 1.0_rk/H * params%gaxx**4 * (-sv_xxaa*neqzp*neqzp/s &
            + sv_xxaa*q(1,:)*q(1,:)*s )&
            + params%gaff*params%gaff*(gam_agff + gam_afgf)/s/H
  !rhs(2,:) = 1.0_rk/H * params%gaxx**4 * (-sv_xxaa*q(2,:)*q(2,:)*s &
  !          + sv_xxaa*neqazp*neqazp/s )&
  !          + params%gaff*params%gaff*(gam_agff + gam_afgf)/s/H
  Ynew = reshape(rhs,(/N/))
end subroutine region3a_eq
