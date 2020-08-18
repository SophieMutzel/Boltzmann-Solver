subroutine allocate_couplings( params )

  implicit none
  type(type_params), intent(inout)                    :: params
  integer(kind=ik)                                    :: i, j, k, nsqrt
  real(kind=rk), dimension(:), allocatable            :: gaxx, kappa

  nsqrt = int(sqrt(real(params%N)))
  allocate(gaxx(nsqrt))
  allocate(kappa(nsqrt))

  call linspace( params%kappa_range(1), params%kappa_range(2), kappa )
  kappa(:) = 10**kappa(:)
  call linspace( params%gaxx_range(1), params%gaxx_range(2), gaxx)
  gaxx(:) = 10**gaxx(:)

  allocate(params%gaff(params%N))
  allocate(params%gaxx(params%N))

  k = 1
  do i=1,nsqrt
    do j=1,nsqrt
      params%gaxx(k) = gaxx(i)
      params%gaff(k) = kappa(j)/gaxx(i)
      k = k+1
    end do
  end do
  deallocate(gaxx, kappa)
end subroutine allocate_couplings
