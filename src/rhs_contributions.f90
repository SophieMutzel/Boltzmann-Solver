subroutine rhs_contributions( N, lz, Y, params, argsint, rhs )
  implicit none
  real(kind=rk), intent(in)                   :: lz
  real(kind=rk), dimension(N), intent(in)     :: Y
  real(kind=rk), dimension(:,:), intent(out)  :: rhs
  type (type_params), intent(in)              :: params
  type (type_argsint), intent(inout)          :: argsint
  integer(kind=ik), intent(in)                :: N
  real(kind=rk), dimension(nrhs,params%N)     :: q
  real(kind=rk)                               :: mx, T, mf, nc, s, H!, hSM, rar, gSM, gHS
  real(kind=rk)                               :: result
  real(kind=rk)                               :: sv_agff, sv_afgf, sv_xxff
  integer(kind=ik)                            :: ier, i, neval!,nd, nr
  real(kind=rk), dimension(params%N)          :: Tprim, sv_aaxx, sv_xxaa, neqzp, neqazp
  real(kind=rk), dimension(params%N)          :: gam_agff, gam_afgf, gam_xxff

  q = reshape(Y,(/nrhs,params%N/))
  mx = params%mx
  T = mx/10**lz
  argsint%T = T
  Tprim(:) = sqrt(params%gaff) * Ta(T,params)
  ! HS interaction
  do i=1,params%N
    ! DM axion interaction
    argsint%g = params%gaxx(i)
    call sigmav( Tprim(i), params, argsint, "aaxx", sv_aaxx(i) )
    call sigmav( Tprim(i), params, argsint, "xxaa", sv_xxaa(i) )
    ! SM axion interaction
    argsint%g = params%gaff(i)
    call gamma_r_new( T, params, argsint, "agff", gam_agff(i) )
    call gamma_r_new( T, params, argsint, "afgf", gam_afgf(i) )
    ! SM DM interaction
    argsint%g = params%gaff(i)*params%gaxx(i)
    call gamma_r_new( T, params, argsint, "xxff", gam_xxff(i) )
  end do
  ! Y_x
  do i=1,params%N
    neqzp(i) = neq(sqrt(params%gaff(i))*Ta(T,params), params%mx, gDM)
    neqazp(i) = neq(sqrt(params%gaff(i))*Ta(T,params), params%ma, ga)
  end do
  s = ent( T, params )
  H = Hub( T, params )

  rhs(1,:) = lz
  rhs(2,:) = sqrt(params%gaff)*Ta(T,params)
  rhs(3,:) = s/H * sv_aaxx * q(2,:) * q(2,:)
  rhs(4,:) = s/H * sv_xxaa * q(1,:) * q(1,:)
  rhs(5,:) = (gam_agff + gam_afgf)/s/H
  rhs(6,:) = gam_xxff/s/H
  rhs(7,:) = sv_xxaa*neqzp*neqzp/s/H

end subroutine rhs_contributions