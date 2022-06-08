subroutine freeze_in_grid(params, argsint)
implicit none
type (type_params), intent(inout)       :: params
type (type_argsint), intent(inout)      :: argsint
real(kind=rk), allocatable              :: gaxxmax(:), gaff(:), Tprime(:,:)
real(kind=rk)                           :: T_start, abserr, rar, Tmx, gaffmin, gaxxmaxax, gaxxmaxSMx,gaffminmin
real(kind=rk)                           :: gaffnew, z, zpdz, atol, conv_eps, dz, eps, rtol, T
real(kind=rk)                           :: epsabs, epsrel, gam_xxff, gam_agff, ffa, sv_aaxx, ma, mx, e,gam_afgf, gam_ahff
integer(kind=ik)                        :: i, ier, neval, nT=51, j
real(kind=rk),dimension(29)             :: mavals
integer(kind=ik)                        :: io_error, istate, nd, k, ierr, itask
TYPE(VODE_OPTS)                         :: OPTIONS
real(kind=rk), dimension(1)             :: Y,Ynew

! lower bound gaff
  allocate(argsint%drhoa(2,size(params%drhoa(1,:))),gaxxmax(nT),gaff(nT))
  allocate(Tprime(2,nT))


  mavals=(/0.00131826, 0.00158489, 0.00251189, 0.00398107, 0.00630957, 0.01, &
        0.0158489, 0.0251189, 0.0398107, 0.0630957, 0.1, 0.158489, 0.251189, &
        0.398107, 0.630957, 1., 1.58489, 2.51189, 3.98107, 6.30957, 10., &
        15.8489, 25.1189, 39.8107, 63.0957, 100., 158.489, 251.189, 398.107/)
  !call linspace(-1.0_rk, 2.7_rk, Tprime(2,:))
  !Tprime(2,:) = 10.0_rk**Tprime(2,:)
  argsint%g = 1.0_rk
  epsabs=1e-30_rk
  epsrel=1e-30_rk
  argsint%drhoa = params%drhoa
  do e=-4.3_rk,-4.0_rk,0.05_rk!e=-2.0_rk,1.0_rk,0.02_rk!i=1,29!!i=23,29!
    ma = 10.0_rk**e!mavals(i)!
    mx = 10.0_rk*ma
    params%ma = ma
    argsint%ma = ma
    params%mx = mx
    argsint%mx = mx
    call gamma_r_new(mx, argsint, "agffth", gam_agff )
    ! inverse decay a->ff
    ffa = gammav(mx, argsint, "affth")
    call gamma_r_new( mx, argsint, "afgfth", gam_afgf )
    call gamma_r_new( mx, argsint, "ahff", gam_ahff )

    gaffmin = log10(sqrt(Hub(mx,0.0_rk)/(gam_agff/neq(mx,ma,ga)+ 2.0_rk*gam_afgf/neq(mx,ma,ga)+ffa)))
    gaffminmin=gaffmin
    gaffnew = gaffmin

    do while (gaffminmin<=gaffnew .and. gaffnew>-10.0_rk)
      gaffnew = gaffnew - 0.05_rk
      params%gaff(1) = 10.0_rk**gaffnew
      Y=0.0_rk
      z = -1.9_rk!log10(0.1_rk)
      atol = 1e-2_rk ! absolute tolerance
      rtol = 1e-2_rk ! relative tolerance

      OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL,H0=params%dz,HMAX=0.001_rk,MXSTEP=1000000)
      itask = 1
      istate = 1
      do while ( z<0.0_rk )!.or. eps > conv_eps )
        zpdz = z + params%dz_plot
        call boltzmann_axion(1, zpdz, Y, Ynew, params, argsint)
        CALL VODE_F90( boltzmann_axion, 1, Y, z, zpdz, itask, &
                    istate, OPTIONS, params, argsint )
        T = mx/10**z
        !write(*,*) z, Y(1)/T/T/T, neq(T,ma,ga)/T/T/T
        if (Y(1)>=neq(mx/10**z,ma,ga)) then
          gaffminmin = gaffnew
        !  write(*,*) gaffminmin
        end if
        !write(*,*) z, Y(1), neq(mx/10**z,mx,gDM)
        !eps = abs((q_tot(2,1,it)-q_tot(2,1,it-1))/q_tot(2,1,it))

        z = zpdz
      end do
    end do



!    call linspace(log10(mx), log10(30.0_rk*mx), Tprime(2,:))
!!    call linspace(log10(mx), log10(mx/10**(-1.9_rk)), Tprime(2,:))
!    Tprime(2,:) = 10.0_rk**Tprime(2,:)
!    gaffminmin=1.0_rk
!    do j=1,nT
!      call gamma_r_new(Tprime(2,j), argsint, "agffth", gam_agff )
!      !write(*,*) gam_agff
!      ! inverse decay a->ff
!      ffa = gammav(Tprime(2,j), argsint, "affth")
!      !write(*,*) ffa
!      call gamma_r_new( Tprime(2,j), argsint, "afgfth", gam_afgf )
!      write(*,*) log10(mx/Tprime(2,j)), Hub(Tprime(2,j),0.0_rk), 10**(2.0_rk*gaffmin)*(gam_agff/neq(Tprime(2,j),ma,ga)+ 2.0_rk*gam_afgf/neq(Tprime(2,j),ma,ga)+ffa)
!      !gaffmin = log10(sqrt(Hub(Tprime(2,j),0.0_rk)/(gam_agff/neq(Tprime(2,j),ma,ga)+ 2.0_rk*gam_afgf/neq(Tprime(2,j),ma,ga)+ffa)))
!      if (gaffmin<gaffminmin) then
!        gaffminmin=gaffmin
!        !write(*,*) gaffmin, Tprime(2,j), mx
!      end if
!    end do

!    call gamma_r_new(mx, argsint, "agffth", gam_agff )
!    ! inverse decay a->ff
!    ffa = gammav(mx, argsint, "affth")
!    call gamma_r_new( mx, argsint, "afgfth", gam_afgf )
!    !gam_afgf = 0.0_rk
!    gaffmin = log10(sqrt(Hub(mx,0.0_rk)/(gam_agff/neq(mx,ma,ga)+ 2.0_rk*gam_afgf/neq(mx,ma,ga)+ffa)))
!rhoeq(ma,ma,ga)
!    call linspace(gaffmin,0.0_rk,gaff)
!    call gamma_r_new( mx, argsint, "xxffth", gam_xxff )
!
!    call sigmav( mx, params, argsint, "aaxx", sv_aaxx )
!
!    do i=1,nT
!      params%gaff(1)=10.0_rk**gaff(i)
!      gaxxmaxSMx = log10(sqrt(Hub(mx,rhoeq(mx,ma,ga))/&
!                  (params%gaff(1)*params%gaff(1)*gam_xxff/neq(mx,mx,ga))))
!!      do j=1,nT
!!        call qags(rhop_over_rho,argsint,T_RH,&
!!                  Tprime(2,j), 1e-5_rk, 1e-5_rk, rar, abserr, neval, ier)
!!        Tprime(1,j) = Ta(Tprime(2,j),params,rar)
!!      end do
!!      call interp_linear(nT, Tprime(1,:),Tprime(2,:),mx, Tmx)
!      ! neq,a(z')
!      gaxxmaxax = log10(sqrt(sqrt(Hub(mx/10.0_rk,rhoeq(mx,mx,gDM)+rhoeq(mx,ma,ga))/&
!                  (sv_aaxx*neq(mx,ma,ga)))))
!      gaxxmax(i) = min(gaxxmaxax,gaxxmaxSMx)
!      !write(*,*) ma, mx, gaff(i), gaxxmax(i)
!    end do
    write(*,*) ma, gaffminmin!gaff(1), maxval(gaxxmax)
  end do

  !deallocate(argsint%drhoa, Tprime, gaxxmax,gaff)

end subroutine freeze_in_grid
