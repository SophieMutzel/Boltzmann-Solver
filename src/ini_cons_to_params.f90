subroutine ini_cons_to_params( params, argsint )
  implicit none

  ! user defined parameter structure
  type (type_params), intent(inout)            :: params
  type (type_argsint), intent(inout)           :: argsint
  real(kind=rk), dimension(:,:), allocatable   :: geff_heff
  character(len=16)                            :: stmx, stma
  integer(kind=ik)                             :: n, i

  write(stmx,'(F12.4)') params%mx
  write(stma,'(F12.4)') params%ma

  !call read_vec(trim(params%ini_direc)//'initial_valuesz01_255.dat', params%initial_values, nrhs)
  !call read_dbl_from_file(trim(params%ini_direc)//'initial_values_Z01_'//trim(adjustl(stmx))//trim(adjustl(stma))//'.dat', params%init_Y)
  call read_matrix_from_file(trim(params%ini_direc)//'geff_heff.txt', geff_heff, 3, .true.)
  !call read_matrix_from_file(trim(params%ini_direc)//'heffHS'//trim(adjustl(stmx))//trim(adjustl(stma))//'.dat', params%heff_HS, 2, .false.)
  !call read_matrix_from_file(trim(params%ini_direc)//'geffHS'//trim(adjustl(stmx))//trim(adjustl(stma))//'.dat', params%geff_HS, 2, .false.)

  n = size(geff_heff,2)
  allocate(params%heff(2,n))
  allocate(params%geff(2,n))
  do i=1,n
    params%heff(1,i)    = geff_heff(1,n+1-i)
    params%heff(2,i)    = geff_heff(3,n+1-i)
    params%geff(1,i)    = geff_heff(1,n+1-i)
    params%geff(2,i)    = geff_heff(2,n+1-i)
  end do

  allocate(params%geff_HS(2,50),params%heff_HS(2,50))
  call heffHS(params, argsint, params%heff_HS)
  ! have to call heffHS first!
  call geffHS(params, argsint, params%geff_HS)
  !params%heff_HS(1,:) = 10**params%heff_HS(1,:)
  !params%geff_HS(1,:) = 10**params%geff_HS(1,:)

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

   !call read_matrix_from_file(trim(params%ini_direc)//'correct/drhoarho59.dat', params%drhoa_rho, 2, .false.)
   call read_matrix_from_file(trim(params%ini_direc)//'rhoprho59.txt', params%rhoa_rho, 2, .false.)
   call read_matrix_from_file(trim(params%ini_direc)//'alphas.dat', params%alpha_s, 2, .false.)

end subroutine ini_cons_to_params
