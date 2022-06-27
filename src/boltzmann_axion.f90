subroutine boltzmann_axion( N, lz, Y, Ynew, params, argsint )
  implicit none
  real(kind=rk), intent(in)                   :: lz
  real(kind=rk), dimension(N), intent(in)     :: Y
  real(kind=rk), dimension(N), intent(inout)  :: Ynew
  type (type_params), intent(in)              :: params
  type (type_argsint), intent(inout)          :: argsint
  integer(kind=ik), intent(in)                :: N
  real(kind=rk), dimension(params%N)          :: q, q10
  real(kind=rk)                               :: mx, ma, T, s, H
  integer(kind=ik)                            :: ier, i, neval,nd,nz
  real(kind=rk), dimension(params%N)          :: neqzp,neqazp
  real(kind=rk), dimension(params%N)          :: gam_agff, gam_afgf, ffa, sv_aaxx, gam_ahff
  real(kind=rk)                               :: rhoeqaTp, rhoeqDMTp

  q = Y
  mx = params%mx
  ma = params%ma
  T = mx/10**lz
  !q10 = 10.0_rk**q
  ! rho,eq,a(T')
  rhoeqaTp = rhoeq(T,ma,ga)
  ! rhoeq,DM(T')
  rhoeqDMTp = rhoeq(T,mx,gDM)
  ! Hubble function
  H = Hub( T, rhoeqaTp+rhoeqDMTp )
  do i=1,params%N
    ! neq,DM(z')
    neqzp(i) = neq(T, mx, gDM)
    ! neq,a(z')
    neqazp(i) = neq(T, ma, ga)
    ! SM axion interaction
    ! argsint%g = params%gaff(i)
    call gamma_r_new( T, argsint, "agffth", gam_agff(i) )
    !call gamma_r_new( T, argsint, "agff", gam_agff(i) )
    gam_agff(i) = gam_agff(i)*params%gaff(i)*params%gaff(i)
    gam_afgf(i) = 0.0_rk

    !nz = size(params%gam_afgf,2)
    !call interp_linear(nz, params%gam_afgf(1,:),params%gam_afgf(2,:),lz, gam_afgf(i))
    !call gamma_r_new( T, argsint, "afgfth", gam_afgf(i) )
    !gam_afgf(i) = gam_afgf(i)*params%gaff(i)*params%gaff(i)!

    call gamma_r_new( T, argsint, "ahff", gam_ahff(i))
    gam_ahff(i) = 0.0_rk!gam_ahff(i)*params%gaff(i)*params%gaff(i)!
    ! inverse decay ff->a
    !ffa(i) = gammav(T, argsint, "affth")*params%gaff(i)*params%gaff(i)
    ffa(i) = 0.0_rk
  end do
    !Ynew(:) =  -l10*3.0_rk*q + l10*(-sv_aaxx*q*q+ &
    !          (gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)*(1.0_rk-q/neqazp))/H

    Ynew(:) =  -l10*3.0_rk*q + l10*((gam_agff + 2.0_rk*gam_afgf + ffa*neqazp + gam_ahff)*(1.0_rk-q/neqazp))/H

              !- 3.0_rk + (-sv_aaxx*neqazp*neqazp/q10+ (gam_agff + 2.0_rk*gam_afgf)*(1.0_rk/q10 - 1.0_rk/neqazp) &
              !+ ffa*(neqazp/q10 - 1.0_rk))/H
    !((-sv_aaxx*neqazp*neqazp)/q10(2,:)+ (gam_agff + 2.0_rk*gam_afgf + ffa*neqazp)*(1.0_rk/q10(2,:) - 1.0_rk/neqazp))/H

end subroutine boltzmann_axion
