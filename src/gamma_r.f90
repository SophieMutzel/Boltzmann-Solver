real(kind=rk) function gamma_r(T,params,sigma,p_f)
  implicit none
  type (type_params), intent(in)      :: params
  real(kind=rk), intent(in)        :: T
  real(kind=rk)                    :: mx, gam
  integer                             :: i

  gam = 0d0

  mx = params%mx
  if (p_f) then
    mfs = params%mf
    ncs = params%ncf
  else
    mfs = params%ma
    ncs = 1d0
  end if
  call qagi( sigma, 4d0*max(mx,mf)*max(mx,mf), 1, epsabs, epsrel, result, abserr, neval, ier )
  kernel = @(s,T,mf,nc) sigma(s,params,mf,nc).*(s-4*mx.^2).*sqrt(s).*besselk(1,sqrt(s)./T);
  gam = zeros(length(mx));

  do i=1,size(mfs)
    gam = gam + T/2/pi^4.*integral(@(s) kernel(s,T,mfs(i),ncs(i)),max(4*mx^2,4*mfs(i)^2),Inf);
  end do
  gamma_r = gam
end function gamma_r
