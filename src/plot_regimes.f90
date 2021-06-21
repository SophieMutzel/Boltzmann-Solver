subroutine plot_regimes()
argsint%g = 1.0_rk
  nsqrt = 100
  allocate( gaxx(nsqrt), kappa(nsqrt) )
  allocate(gaxx_tot(nsqrt*nsqrt),gaff_tot(nsqrt*nsqrt))
  call linspace( -18.0_rk, -1.0_rk, kappa )
  kappa(:) = 10**kappa(:)
  call linspace( -12.0_rk, -1.0_rk, gaxx)
  gaxx(:) = 10**gaxx(:)

  k = 1
  do i=1,nsqrt
    do j=1,nsqrt
      gaxx_tot(k) = gaxx(i)
      gaff_tot(k) = kappa(j)/gaxx(i)
      k = k+1
    end do
  end do
open (unit=96, file="temp/regimes2.txt", status='old', action='write', position='append', iostat=io_error)
if (io_error==0) then
do i=1,size(gaxx_tot)
  params%gaff(1) = gaff_tot(i)
  params%gaxx(1) = gaxx_tot(i)
  call choose_regime(params, argsint)
    select case(params%regime)
    case("freeze-in")
      write(96,*) params%gaxx(1)*params%gaff(1), params%gaxx(1), 1
    case("seq-freeze-in")
      write(96,*) params%gaxx(1)*params%gaff(1), params%gaxx(1), 2
    case("reannihilation")
      write(96,*) params%gaxx(1)*params%gaff(1), params%gaxx(1), 3
    case("freeze-out")
      write(96,*) params%gaxx(1)*params%gaff(1), params%gaxx(1), 4
    end select
end do
  else
    WRITE(*,*) "error"
  end if
  close(96)

end subroutine plot_regimes
