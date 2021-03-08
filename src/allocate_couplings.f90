subroutine allocate_couplings( params )

  implicit none
  type(type_params), intent(inout)            :: params
  integer(kind=ik)                            :: i, j, k, nsqrt, ierr, nmin, nextra
  integer(kind=ik), dimension(:), allocatable :: sendcounts, displs
  real(kind=rk), dimension(:), allocatable    :: gaxx, kappa, gaxx_tot, gaff_tot

  allocate( gaxx_tot(params%N_tot), gaff_tot(params%N_tot) )
  allocate( sendcounts(0:params%nprocs-1), displs(0:params%nprocs-1) )
  if ( params%rank == 0 ) then
    nsqrt = int(sqrt(real(params%N_tot)))
    allocate( gaxx(nsqrt), kappa(nsqrt) )

    call linspace( params%kappa_range(1), params%kappa_range(2), kappa )
    kappa(:) = 10**kappa(:)
    call linspace( params%gaxx_range(1), params%gaxx_range(2), gaxx)
    gaxx(:) = 10**gaxx(:)

    k = 1
    do i=1,nsqrt
      do j=1,nsqrt
        gaxx_tot(k) = gaxx(i)
        gaff_tot(k) = kappa(j)/gaxx(i)
        k = k+1
      end do
    end do
  end if

  nmin = params%N_tot/params%nprocs
  nextra = mod(params%N_tot,params%nprocs)
  k = 0
  do i = 0, params%nprocs-1
     if (i < nextra) then
        sendcounts(i) = nmin + 1
     else
        sendcounts(i) = nmin
     end if
     displs(i) = k
     k = k + sendcounts(i)
   end do
   params%N  = sendcounts(params%rank)


   allocate(params%gaff(params%N))
   allocate(params%gaxx(params%N))
   ! distribute couplings among different ranks
   call MPI_Scatterv( gaxx_tot, sendcounts, displs, MPI_DOUBLE, &
       params%gaxx, params%N, MPI_DOUBLE, &
       0, params%BOLTZMANN_COMM, ierr)
   call MPI_Scatterv( gaff_tot, sendcounts, displs, MPI_DOUBLE, &
       params%gaff, params%N, MPI_DOUBLE, &
       0, params%BOLTZMANN_COMM, ierr)

  params%gaxx=10**params%gaxx_range(1)!0.01_rk!0.005_rk!0.00501187_rk
  params%gaff=10**(params%kappa_range(1)-params%gaxx_range(1))!1.0e-12_rk!6.30957e-12_rk!1.99526e-13_rk!
  if ( params%rank == 0 ) deallocate( gaxx, kappa )

  deallocate( gaxx_tot, gaff_tot, sendcounts, displs )
end subroutine allocate_couplings
