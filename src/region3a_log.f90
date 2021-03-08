subroutine region3a_log( N, lz, Y, Ynew, params, argsint )
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
  real(kind=rk), dimension(params%N)                :: Tprim, sv_aaxx, sv_xxaa

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  T = mx/10**lz
  !T= mx/lz
  argsint%T = T
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
  s = ent( T, params )
  H = Hub( T, params )
  ! Y_x
  rhs(1,:) = s/H *params%gaxx**4 * (sv_aaxx*q(2,:)*q(2,:)&
            - sv_xxaa*q(1,:)*q(1,:))!+&
            !params%gaff*params%gaff*params%gaxx*params%gaxx*&
            !gam_xxff/s/H
  ! Y_a
  rhs(2,:) = params%gaxx**4 * (-sv_aaxx*q(2,:)*q(2,:) &
            + sv_xxaa*q(1,:)*q(1,:))*s/H !&
            !+ params%gaff*params%gaff*(gam_agff + gam_afgf)/s/H
  Ynew = reshape(rhs,(/N/))
end subroutine region3a_log
