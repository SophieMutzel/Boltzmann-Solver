subroutine ini_cons_to_params( params )
  implicit none

  ! user defined parameter structure
  type (type_params), intent(inout)            :: params
  real(kind=rk), dimension(:,:), allocatable   :: geff_heff
  character(len=16)                            :: stmx, stma
  integer(kind=ik)                             :: n, i

  write(stmx,'(F12.4)') params%mx
  write(stma,'(F12.4)') params%ma

  call read_dbl_from_file(trim(params%ini_direc)//'initial_values_01_'//trim(adjustl(stmx))//trim(adjustl(stma))//'.dat', params%init_rhoprime)
  call read_dbl_from_file(trim(params%ini_direc)//'initial_values_Z01_'//trim(adjustl(stmx))//trim(adjustl(stma))//'.dat', params%init_Y)
  call read_matrix_from_file(trim(params%ini_direc)//'geff_heff.txt', geff_heff, 3,.true.)
  call read_matrix_from_file(trim(params%ini_direc)//'heffHS'//trim(adjustl(stmx))//trim(adjustl(stma))//'.dat', params%heff_HS, 2,.false.)
  call read_matrix_from_file(trim(params%ini_direc)//'geffHS'//trim(adjustl(stmx))//trim(adjustl(stma))//'.dat', params%geff_HS, 2,.false.)

  n = size(geff_heff,2)
  allocate(params%heff(2,n))
  allocate(params%geff(2,n))
  do i=1,n
    params%heff(1,i)    = geff_heff(1,n+1-i)
    params%heff(2,i)    = geff_heff(3,n+1-i)
    params%geff(1,i)    = geff_heff(1,n+1-i)
    params%geff(2,i)    = geff_heff(2,n+1-i)
  end do

  params%heff_HS(1,:) = 10**params%heff_HS(1,:)
  params%geff_HS(1,:) = 10**params%geff_HS(1,:)

end subroutine ini_cons_to_params
