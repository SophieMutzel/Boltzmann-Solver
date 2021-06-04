subroutine ini_cons_to_params( params, argsint )
  implicit none

  ! user defined parameter structure
  type (type_params), intent(inout)            :: params
  type (type_argsint), intent(inout)           :: argsint
  real(kind=rk), dimension(:,:), allocatable   :: fgfa
  character(len=16)                            :: stmx, stma
  integer(kind=ik)                             :: n, i, nc
  real(kind=rk)                                :: xmin,xmax,ymin,ymax

  write(stmx,'(F12.4)') params%mx
  write(stma,'(F12.4)') params%ma

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

   !call read_matrix_from_file(trim(params%ini_direc)//'rhoprho59.txt', params%rhoa_rho, 2, .false.)
   call read_matrix_from_file(trim(params%ini_direc)//'drhop09thlogT.txt', params%drhoa, 2, .false.)
   call read_matrix_from_file(trim(params%ini_direc)//'drhopfgfa09th.txt', fgfa, 2, .false.)
   params%drhoa(2,:) = params%drhoa(2,:) + fgfa(2,:)
   deallocate(fgfa)
   !call read_matrix_from_file(trim(params%ini_direc)//'drhop09.txt', params%drhoa, 2, .false.)

end subroutine ini_cons_to_params
