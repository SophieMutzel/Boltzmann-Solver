subroutine rhs_boltzmann( q, params, argsint, z, rhs )

  implicit none
  real(kind=rk), intent(in)                         :: z
  real(kind=rk), dimension(:,:), intent(in)         :: q
  real(kind=rk), dimension(:,:), intent(inout)      :: rhs
  type (type_params), intent(in)                    :: params
  type (type_argsint), intent(inout)                :: argsint
  real(kind=rk)                                     :: mx, g1, g2, T, gHS, mf, nc, hSM
  real(kind=rk)                                     :: result, gamma_rcon
  real(kind=rk)                                     :: epsabs, epsrel, abserr
  integer(kind=ik)                                  :: ier, i, nd, neval
  real(kind=rk), dimension(params%N)                :: integral, Tprim, gamma_rprim, Yeqprim

  epsabs = 1e-5_rk
  epsrel = 1e-5_rk
  mx = params%mx
  g1 = 4.0_rk
  g2 = 4.0_rk

  T = mx/z
  argsint%T = T
  nd = size(params%heff_HS,2)
  call interp_linear(nd, params%heff_HS(1,:),params%heff_HS(2,:),T, gHS)
  Tprim(:) = ( heff( T, params )/gHS * q(2,:) )**(0.25_rk)*T;

  integral = 0.0_rk

  do i=1,size(params%mf)
    argsint%mf = params%mf(i)
    argsint%nc = params%ncf(i)
    call qagi( kernel, argsint, 4.0_rk*max(mx*mx,params%mf(i)*params%mf(i)), 1,&
                epsabs, epsrel, result, abserr, neval, ier )
    integral = integral + params%gaff*params%gaff*params%gaxx*params%gaxx * result
  end do

  do i=1,params%N
    call gamma_r( Tprim(i), params, argsint, .false., gamma_rprim(i))
    Yeqprim(i) = Yeq( Tprim(i), params )
  end do

  call gamma_r( T, params, argsint, .true., gamma_rcon)
  
  rhs(1,:) = 1.0_rk/( Hub( T, params )* z *ent( T, params ) )* &
          ( params%gaff*params%gaff*params%gaxx*params%gaxx * gamma_rcon* &
          ( 1.0_rk - ( q(1,:) / Yeq(T,params) ) * ( q(1,:) / Yeq(T,params) ) ) &
          + params%gaxx**4 * gamma_rprim * &
          ( 1.0_rk - ( q(1,:) / Yeqprim ) * ( q(1,:) / Yeqprim )) )

  call heffSM( T, params, hSM )
  rhs(2,:) = z*z*z* 15.0_rk / ( 16.0_rk * Hub( T, params )* pi**6 *hSM *mx**4 )*g1*g2*integral

end subroutine rhs_boltzmann
