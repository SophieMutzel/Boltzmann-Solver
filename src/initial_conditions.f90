! intitial conditions, different depending on regime we are in
subroutine initial_conditions( params, q, q_tot, rhs, argsint )
  implicit none
  type (type_params), intent(inout)                              :: params
  type (type_argsint), intent(inout)                             :: argsint
  real(kind=rk), dimension(:,:), allocatable, intent(inout)      :: q
  real(kind=rk), dimension(:,:,:), allocatable, intent(inout)    :: q_tot, rhs
  integer(kind=ik)                                               :: L, i, nr, ier, neval
  real(kind=rk)                                                  :: T_start, s, rar, rho, Tprime, gam_agff, ffa, aaxx
  real(kind=rk)                                                  :: epsabs, epsrel, abserr, Ya, rhoovern_a, ma,j,k,T

  allocate(q(nrhs,params%N))
  select case (params%regime)
    case ("reannihilation")
      allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
      epsabs=1e-30_rk
      epsrel=1e-30_rk
      params%z_start = log10(params%mx/T_RH)+0.1_rk!-2.5_rk
      Tprime=1000000.0_rk
      argsint%drhoa = params%drhoa
      do while (Tprime>10.0_rk*params%mx)
        params%z_start = params%z_start+0.2_rk
        T_start = params%mx/10**params%z_start
        write(*,*) T_start
        call qags(rhop_over_rho,argsint,T_RH,&
              T_start, epsabs, epsrel, rar, abserr, neval, ier)
        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
        Tprime = Ta(T_start,params,rar)
      end do
      if (Tprime>T_start)  then
        write(*,*) "T= ", T_start, "T'= ", Tprime, params%z_start
        write(*,*) "WARNING: T'>T in reannihilation regime! Setting T' to T. "
        Tprime = T_start
        call abort_it("T'>T in reannihilation regime!")
      end if
      do i=1,params%N
        ! Y_chi=Yeq(T')
        q(1,i) = neq(Tprime, params%mx, gDM)
        ! Y_a=Yeq,a(T')
        q(2,i) = neq(Tprime, params%ma, ga)
        q(3,i) = Tprime
        write(*,*) q, params%z_start
      end do
      params%z_max =  3.0_rk
      deallocate(argsint%drhoa)
    case ("seq-freeze-in")
      params%z_start = -1.9_rk
      params%z_max = 1.0_rk
      T_start = params%mx/10**params%z_start
      Tprime = T_start
      ! Y_chi=0
      q(1:2,:) = 1e-30_rk
      do i=1,params%N
        q(3,i) = Tprime
      end do
    case ("freeze-in")
      params%z_start = -1.9_rk
      params%z_max = 1.0_rk
      T_start = params%mx/10**params%z_start
      Tprime = T_start
      ! Y_chi=0
      q(1:2,:) = 1e-30_rk
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
  L = int((params%z_max-params%z_start)/params%dz_plot)+50

  allocate(q_tot(nrhs+5,params%N,L+5))
  allocate(rhs(nrhs+6,params%N,L+5))

  s = T_start*T_start*T_start!ent(T_start, params)
  q_tot(1,:,1) = params%z_start
  q_tot(2:nrhs,:,1) = q(1:2,:)/s
  q_tot(nrhs+1,:,1) = neq(T_start, params%mx, gDM)/s
  q_tot(nrhs+2,:,1) = neq(Tprime, params%mx, gDM)/s!q(1,:)/s
  q_tot(nrhs+3,:,1) = neq(Tprime, params%ma, ga)/s!q(2,:)/s
  q_tot(nrhs+4,:,1) = Tprime
  q_tot(2:nrhs,:,1) = log10(q(1:2,:)/s)
  q(1:2,:) = log10(q(1:2,:))


end subroutine initial_conditions
