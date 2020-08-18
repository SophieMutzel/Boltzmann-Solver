
module module_params
  use module_precision
  use module_read_write
  use module_utils

  implicit none

  ! global user defined data structure for time independent variables
  type type_params

    ! masses
    real(kind=rk)                                :: mx = 0.25, ma = 0.05
    real(kind=rk), dimension(9)                  :: mf = (/ 0.000511, 0.1057, 1.777, 0.0046, 0.096, 4.18, 0.0022, 1.28, 173.0 /)
    real(kind=rk), dimension(9)                  :: ncf = (/ 1., 1., 1., 3., 3., 3., 3., 3., 3. /)

    ! couplings
    integer(kind=ik)                             :: N = 16
    real(kind=rk), dimension(:), allocatable     :: gaff, gaxx
    real(kind=rk), dimension(2)                  :: kappa_range = (/ -10., -8. /), gaxx_range = (/ -10., -8. /)

    ! constants
    real(kind=rk)                                :: Mpl = 1.22e19, g = 10.75

    ! start, end of simulation
    real(kind=rk)                                :: z_max = 100.0, z_start = 0.01, dz = 0.00001
    integer(kind=ik)                             :: nt = 100

    ! correct relic density?
    real(kind=rk)                                :: DM_low = 3.8e-10, DM_up = 4.3e-10

    ! initial conditions
    real(kind=rk)                                :: init_Y, init_rhoprime
    character(len=75)                            :: ini_direc="/Users/sophiemutzel 1/Documents/Wolfram Mathematica/Axion_DM/RelicDensity/"

    ! effective degrees of freedom
    real(kind=rk), dimension(:,:), allocatable   :: heff, geff, heff_HS, geff_HS

end type type_params

contains
    include "allocate_couplings.f90"
    include "ini_cons_to_params.f90"
    include "initial_conditions.f90"

end module module_params
