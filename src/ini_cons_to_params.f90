subroutine ini_cons_to_params( params, argsint )
  implicit none

  ! user defined parameter structure
  type (type_params), intent(inout)            :: params
  type (type_argsint), intent(inout)           :: argsint
  real(kind=rk), dimension(:,:), allocatable   :: geff_heff, points
  character(len=16)                            :: stmx, stma
  integer(kind=ik)                             :: n, i, nc
  real(kind=rk)                                :: xmin,xmax,ymin,ymax

  write(stmx,'(F12.4)') params%mx
  write(stma,'(F12.4)') params%ma

  !call read_matrix_from_file(trim(params%ini_direc)//'geff_heff.txt', geff_heff, 3, .true.)
!  call read_matrix_from_file(trim(params%ini_direc)//'geff_heff_new.txt', geff_heff, 5, .true.)
!
!  n = size(geff_heff,2)
!  allocate(params%heff(2,n))
!  allocate(params%geff(2,n))
!  params%geff(1,:) = geff_heff(1,:)
!  params%heff(1,:) = geff_heff(1,:)
!
!  params%geff(2,:) = geff_heff(2,:)!+geff_heff(3,:)
!  params%heff(2,:) = geff_heff(4,:)!+geff_heff(5,:)

  allocate(params%geff_HS(2,100),params%heff_HS(2,100))!,params%B(500))!,params%C(500),params%D(500))
  call heffHS(params, argsint, params%heff_HS)
  ! have to call heffHS first!
  call geffHS(params, argsint, params%geff_HS)

  ! for bezier fit to geff
!  nc = size(params%geff(1,:))
!  n = size(params%geff(1,1:nc:3))
!  allocate(params%A(n-1,2),params%B(n-1,2),points(n,2))
!  xmin = minval(params%geff(1,1:nc:3))
!  xmax = maxval(params%geff(1,1:nc:3))
!  ymin = minval(params%geff(2,1:nc:3))
!  ymax = maxval(params%geff(2,1:nc:3))
!  points(:,1) = (params%geff(1,1:nc:3)-xmin)/(xmax-xmin)
!  points(:,2) = (params%geff(2,1:nc:3)-ymin)/(ymax-ymin)
!  call get_bezier_coeff(points,n-1,params%A,params%B)
!  deallocate(points)

 ! loop factor (sum_f mf^2 ncf qf^2 C0(0,0,ma^2,mf^2, mf^2,mf^2))^2
!  select case (params%ma)
!
!    case (5.0_rk)
!      params%C0sq = 5.83845_rk
!    case (0.005_rk)
!      params%C0sq = 7.80685_rk
!    case (0.00005_rk)
!      params%C0sq = 16.3943_rk
!    case default
!      write(*,*) "Loop coupling not known for ma=", params%ma
!   end select

   call read_matrix_from_file(trim(params%ini_direc)//'rhoprho59.txt', params%rhoa_rho, 2, .false.)
   call read_matrix_from_file(trim(params%ini_direc)//'drhop59.txt', params%drhoa, 2, .false.)
   call read_matrix_from_file(trim(params%ini_direc)//'alphas.dat', params%alpha_s, 2, .false.)

end subroutine ini_cons_to_params
