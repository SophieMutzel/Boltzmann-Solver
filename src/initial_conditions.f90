subroutine initial_conditions( params, q, q_tot, rhs, argsint )
  implicit none
  type (type_params), intent(inout)                              :: params
  type (type_argsint), intent(inout)                             :: argsint
  real(kind=rk), dimension(:,:), allocatable, intent(inout)      :: q
  real(kind=rk), dimension(:,:,:), allocatable, intent(inout)    :: q_tot, rhs
  integer(kind=ik)                                               :: L, i, nr, ier, neval
  real(kind=rk)                                                  :: T_start, s, rar, rho,Tprime
  real(kind=rk)                                                  :: epsabs, epsrel, abserr

  allocate(q(nrhs,params%N))

  select case (params%regime)
    case ("reannihilation")
      params%z_start = -1.0_rk
      params%z_max = 2.0_rk
      T_start = params%mx/10**params%z_start
      epsabs=1e-30_rk
      epsrel=1e-30_rk
      allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
      argsint%drhoa = params%drhoa
      call qags(rhop_over_rho,argsint,T_RH,&
              T_start, epsabs, epsrel, rar, abserr, neval, ier)
      if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
      deallocate(argsint%drhoa)
      Tprime = Ta(T_start,params,rar)
      do i=1,params%N
        ! Y_chi=Yeq(T')
        q(1,i) = neq(Tprime, params%mx, gDM)
        ! Y_a=Yeq,a(T')
        q(2,i) = neq(Tprime, params%ma, ga)
        q(3,i) = Tprime
      end do
    case ("seq-freeze-in")
      !Tprime = 0.0_rk
      params%z_start = -1.8_rk
      params%z_max = 0.5_rk
      T_start = params%mx/10**params%z_start
      epsabs=1e-30_rk
      epsrel=1e-30_rk
      allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
      argsint%drhoa = params%drhoa
      call qags(rhop_over_rho,argsint,T_RH,&
              T_start, epsabs, epsrel, rar, abserr, neval, ier)
      if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
      deallocate(argsint%drhoa)
      Tprime = Ta(T_start,params,rar)
      do i=1,params%N
        ! Y_chi=Yeq(T')
        q(1,i) = 0.0_rk
        ! Y_a=Yeq,a(T')
        q(2,i) = neq(Tprime, params%ma, ga)
        q(3,i) = Tprime
      end do
    case ("freeze-in")
      params%z_start = -1.8_rk
      params%z_max = 0.5_rk
      T_start = params%mx/10**params%z_start
      Tprime = T_start
      ! Y_chi=0
      q(1:2,:) = 0.0_rk
      do i=1,params%N
        ! Y_a=Yeq,a(T')
        !q(2,i) = neq(Tprime, params%ma, ga)
        q(3,i) = Tprime
      end do
    case ("freeze-out")
      params%z_start = 0.2_rk
      params%z_max = 2.0_rk
      T_start = params%mx/10**params%z_start
      Tprime = T_start
      do i=1,params%N
        ! Y_chi=Yeq(T)
        q(1,i) = neq(Tprime, params%mx, gDM)
        ! Y_a=Yeq,a(T)
        q(2,i) = neq(Tprime, params%ma, ga)
        q(3,i) = Tprime
      end do
    case default
      write(*,*) "Error! Regime ", params%regime, " not (yet) implemented!"
  end select

  L = int((params%z_max-params%z_start)/params%dz_plot)
  allocate(q_tot(nrhs+5,params%N,L+3))
  allocate(rhs(nrhs+5,params%N,L+3))

  s = ent(T_start, params)

  q_tot(1,:,1) = params%z_start
  q_tot(2:nrhs,:,1) = q(1:2,:)/s
  q_tot(nrhs+1,:,1) = neq(T_start, params%mx, gDM)/s
  q_tot(nrhs+2,:,1) = q(1,:)/s
  q_tot(nrhs+3,:,1) = q(2,:)/s
  q_tot(nrhs+4,:,1) = Tprime

end subroutine initial_conditions
