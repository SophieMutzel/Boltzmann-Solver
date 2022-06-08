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

    ! regime
    character(len=14)                            :: regime

    ! couplings
    integer(kind=ik)                             :: N_tot = 1, N
    real(kind=rk), dimension(:), allocatable     :: gaff, gaxx
    real(kind=rk), dimension(2)                  :: kappa_range = (/ -13.5, -13.5 /), gaxx_range = (/ -2.18251, -2.18251 /)

    ! start, end of simulation
    real(kind=rk)                                :: z_max = 2.0_rk, z_start = -1.8, dz = 0.000000001_rk, dz_plot=0.1_rk!z_max = 100.0, z_start = 0.01, dz = 0.01
    integer(kind=ik)                             :: nt = 100000

    ! initial conditions
    real(kind=rk), dimension(3)                  :: initial_values
    character(len=75)                            :: ini_direc="/Users/sophiemutzel/Documents/Wolfram Mathematica/Axion_DM/RelicDensity/"

    character(len=80)                            :: file
    ! effective degrees of freedom
    real(kind=rk), dimension(:,:), allocatable   :: heff, geff, heff_HS, geff_HS

    ! in Mathematica calculated rhs for rhoa/rho
    real(kind=rk), dimension(:,:), allocatable   :: drhoa

    ! in Mathematica calculated rhoa/rho
    real(kind=rk), dimension(:,:), allocatable   :: rhoa_rho

    ! precalculate gamma(af->gf) to accelerate computation
    real(kind=rk), dimension(:,:), allocatable   :: gam_afgf

    logical                                      :: grid = .false.

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
