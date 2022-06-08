subroutine aftogf(params, argsint)
  implicit none
  type (type_params), intent(inout)        :: params
  type (type_argsint), intent(inout)       :: argsint
  real(kind=rk)                            :: T
  integer(kind=ik)                         :: j, nz=200

  allocate(params%gam_afgf(2,nz))
  !call linspace(params%z_start,params%z_max,params%gam_afgf(1,:))
  !do j=1,nz
  !  T = params%mx/10**params%gam_afgf(1,j)
  !  call gamma_r_new( T, argsint, "afgfth", params%gam_afgf(2,j) )
  !end do

  call read_matrix_from_file("/Users/sophiemutzel/Desktop/Uni/Masterarbeit/Boltzmann_solver/afgfma1mx10new.txt", params%gam_afgf, 2, .false.)

end subroutine aftogf
