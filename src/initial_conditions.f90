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
!  allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
!  argsint%g = 1.0_rk
!  epsabs=1e-30_rk
!  epsrel=1e-30_rk
!  argsint%drhoa = params%drhoa
!  do j=-15,-9
!    do k=-3,1,0.5
!      params%gaff(1)=10.0_rk**j
!      params%ma=10.0_rk**k
!      params%mx=10.0_rk*params%ma
!      argsint%ma=params%ma
!      argsint%mx=params%mx
!      params%z_start = k-1.0_rk
!      Tprime=1000000.0_rk
!      do while (Tprime>params%mx)
!        params%z_start = params%z_start+0.2_rk
!        T_start = params%mx/10**params%z_start
!        call qags(rhop_over_rho,argsint,T_RH,&
!              T_start, epsabs, epsrel, rar, abserr, neval, ier)
!        !if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
!        Tprime = Ta(T_start,params,rar)
!      end do
!      !write(*,*) Tprime, params%mx
!      call sigmav( Tprime, params, argsint, "aaxx", aaxx )
!      write(*,*) params%ma, params%gaff(1), &
!      sqrt(sqrt(Hub(T_start,rhoeq(Tprime,params%mx,gDM)+rhoeq(Tprime,params%ma,ga))/aaxx/neq(Tprime,params%ma,ga)))
!    end do
!  end do
!  stop
!  argsint%g = 1.0_rk
!  allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
!  argsint%drhoa = params%drhoa
!  write(*,*) argsint%drhoa(1,:)
!  do j=-3,3,0.1
!    !ma=10.0_rk**j
!    T=params%mx/10.0_rk**j
!    argsint%ma=10.0_rk
!    argsint%mx=100.0_rk
!    epsabs=1e-30_rk
!    epsrel=1e-30_rk
!    call qags(rhop_over_rho,argsint,T_RH,&
!          T, epsabs, epsrel, rar, abserr, neval, ier)
!    write(*,*) j, rar!sqrt(1.0_rk/rar)
!    !sqrt(0.5_rk/rar/rho_SM(0.001_rk)*7.0_rk/8.0_rk*ggamma*pi*pi/30.0_rk*1e-12_rk)
!  end do
!  stop

  select case (params%regime)
    case ("reannihilation")
      allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
      epsabs=1e-30_rk
      epsrel=1e-30_rk
      params%z_start = -2.5_rk
      Tprime=1000000.0_rk
      argsint%drhoa = params%drhoa
      do while (Tprime>params%mx)
        params%z_start = params%z_start+0.2_rk
        T_start = params%mx/10**params%z_start
        call qags(rhop_over_rho,argsint,T_RH,&
              T_start, epsabs, epsrel, rar, abserr, neval, ier)
        if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition"
        Tprime = Ta(T_start,params,rar)
      end do
!      argsint%ra_ini = 0.0_rk
!      argsint%g = params%gaff(1)
!      call qags(n_axion_seq, argsint,T_start,&
!              1000.0_rk, epsabs, epsrel, Ya, abserr, neval, ier)
!      rhoovern_a = Ya/ent(T_start,params)
!      argsint%g = 1.0_rk
!      call gamma_r_new( params%ma, argsint, "agffth", gam_agff )
!      gam_agff = gam_agff*params%gaff(1)*params%gaff(1)
!      ! inverse decay a->ff
!      ffa = gammav(T_start, argsint, "affth")*params%gaff(1)*params%gaff(1)
!      write(*,*) (gam_agff/neq(T_start,params%ma,ga)**2*rhoovern_a+ffa)/Hub(T_start,0.0_rk)
      if (Tprime>T_start)  then
        write(*,*) "T= ", T_start, "T'= ", Tprime, params%z_start
        write(*,*) "WARNING: T'>T in reannihilation regime! Setting T' to T. "
        Tprime = T_start
        !call abort_it("T'>T in reannihilation regime!")
      end if
      do i=1,params%N
        ! Y_chi=Yeq(T')
        q(1,i) = neq(Tprime, params%mx, gDM)
        ! Y_a=Yeq,a(T')
        q(2,i) = neq(Tprime, params%ma, ga)
        q(3,i) = Tprime
      end do
      params%z_max =  3.0_rk
      deallocate(argsint%drhoa)
    case ("seq-freeze-in")
!      params%z_start = -1.8_rk
!      params%z_max = 0.5_rk
!      T_start = params%mx/10**params%z_start
!      epsabs=1e-30_rk
!      epsrel=1e-30_rk
!      allocate(argsint%drhoa(2,size(params%drhoa(1,:))))
!      argsint%drhoa = params%drhoa
!      argsint%g = params%gaff(1)
!      call qags(rhop_over_rho,argsint,T_RH,&
!              T_start, epsabs, epsrel, rar, abserr, neval, ier)
!      if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition T'"
!      deallocate(argsint%drhoa)
!      argsint%ra_ini = rar*rho_SM(T_start)!*params%gaff(1)*params%gaff(1)
!      call qags(n_axion_seq, argsint,T_start,&
!              T_RH, epsabs, epsrel, Ya, abserr, neval, ier)
!      if (ier > 0) write(*,*) "Integral did not converge ier=", ier, " err=", abserr, "for initial condition na"
!      rhoovern_a = rar*rho_SM(T_start)/Ya/ent(T_start,params)
!      Tprime = Ta_seq_fi(T_start,params,rhoovern_a)
!      do i=1,params%N
!        ! Y_chi=Yeq(T')
!        q(1,i) = 0.0_rk
!        ! Y_a=Yeq,a(T')
!        q(2,i) = Ya*ent(T_start,params)!*params%gaff(i)*params%gaff(i)
!        q(3,i) = Tprime
!      end do
!      argsint%g = 1.0_rk
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
  allocate(rhs(nrhs+5,params%N,L+5))

  s = T_start*T_start*T_start!ent(T_start, params)
  q_tot(1,:,1) = params%z_start
  q_tot(2:nrhs,:,1) = q(1:2,:)/s
  q_tot(nrhs+1,:,1) = neq(T_start, params%mx, gDM)/s
  q_tot(nrhs+2,:,1) = neq(Tprime, params%mx, gDM)/s!q(1,:)/s
  q_tot(nrhs+3,:,1) = neq(Tprime, params%ma, ga)/s!q(2,:)/s
  q_tot(nrhs+4,:,1) = Tprime
!  if (params%regime == "reannihilation" .or. params%regime == "freeze-out") then
    q_tot(2:nrhs,:,1) = log10(q(1:2,:)/s)
    q(1:2,:) = log10(q(1:2,:))

!  else
!    q_tot(2,:,1) = log10(1e-40_rk/s)
!    q_tot(3,:,1) = log10(q(2,:)/s)
!  end if


end subroutine initial_conditions
