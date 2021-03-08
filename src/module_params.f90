module module_params

  use mpi
  use module_precision
  use module_read_write
  use module_utils

  implicit none

  ! global user defined data structure for time independent variables
  type, public :: type_params

    ! masses
    real(kind=rk)                                :: mx = 25.0_rk, ma = 5.0_rk

    ! loop coefficient
    real(kind=rk)                                :: C0sq

    ! couplings
    integer(kind=ik)                             :: N_tot = 1, N
    real(kind=rk), dimension(:), allocatable     :: gaff, gaxx
    real(kind=rk), dimension(2)                  :: kappa_range = (/ -13.5, -13.5 /), gaxx_range = (/ -2.18251, -2.18251 /)

    ! constants
    real(kind=rk)                                :: Mpl = 1.22e19, g = 10.75

    ! start, end of simulation
    real(kind=rk)                                :: z_max = 2.0_rk, z_start = -1.8, dz = 0.000000001_rk, dz_plot=0.1_rk!z_max = 100.0, z_start = 0.01, dz = 0.01
    integer(kind=ik)                             :: nt = 100000

    ! correct relic density?
    real(kind=rk)                                :: DM_low = 1.9e-10, DM_up = 2.2e-10

    ! initial conditions
    real(kind=rk), dimension(3)                  :: initial_values
    character(len=75)                            :: ini_direc="/Users/sophiemutzel 1/Documents/Wolfram Mathematica/Axion_DM/RelicDensity/"

    character(len=80)                            :: file
    ! effective degrees of freedom
    real(kind=rk), dimension(:,:), allocatable   :: heff, geff, heff_HS, geff_HS

    ! in Mathematica calculated rhs for rhoa/rho
    real(kind=rk), dimension(:,:), allocatable   :: drhoa_rho

    ! in Mathematica calculated rhoa/rho
    real(kind=rk), dimension(:,:), allocatable   :: rhoa_rho

    ! running of strong coupling alpha_s
    real(kind=rk), dimension(:,:), allocatable   :: alpha_s

    ! process rank
    integer(kind=ik)                            :: rank = -1
    ! number of processes
    integer(kind=ik)                            :: nprocs = -1
    ! WABBIT communicator
    integer(kind=ik)                            :: BOLTZMANN_COMM

end type type_params

contains
    include "allocate_couplings.f90"
    include "get_params.f90"

end module module_params
