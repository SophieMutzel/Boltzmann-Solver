! constants
module module_precision
  use mpi
  implicit none
  integer, parameter, public         :: rk=selected_real_kind(8)
  integer, parameter, public         :: ik=selected_int_kind(8)

  real(kind=rk), parameter, public   :: pi  = 4.0 * atan(1.0)
  real(kind=rk), parameter, public   :: l10 = log(10.0_rk)
  integer(kind=ik)                   :: BOLTZMANN_COMM
  ! dofs of ALP, photon, gluon, DM
  real(kind=rk),parameter,public     :: ga=1.0_rk, ggamma=2.0_rk, gg=16.0_rk, gDM=2.0_rk
  ! SM fermion masses, in GeV
  real(kind=rk), dimension(9)        :: mf = (/ 0.000511, 0.1057, 1.777, 0.0046, 0.096, 4.18, 0.0022, 1.28, 173.0 /)
  ! SM colour dof
  real(kind=rk), dimension(9)        :: ncf = (/ 1., 1., 1., 3., 3., 3., 3., 3., 3. /)
  ! dof of SM fermions
  real(kind=rk), dimension(9)        :: gf = (/2., 2., 2., 6., 6., 6., 6., 6., 6. /)
  ! electric charges of SM fermions
  real(kind=rk), dimension(9)        :: qf = (/ 1., 1., 1., -1./3.,-1./3.,-1./3.,2./3.,2./3.,2./3. /)
  ! reheating temperature
  real(kind=rk), parameter, public   :: T_RH=2000.0_rk
  ! fine structure constant
  real(kind=rk), parameter, public   :: alpha_QED = 1.0_rk/137.0_rk
  ! Planck mass
  real(kind=rk), parameter, public   :: Mpl = 1.22e19_rk
  ! observed x density today*m_x (notice that in our model mDM=2mx)
  real(kind=rk), parameter, public   :: Yxmxobs = 2.045e-10_rk
  ! Below 600 MeV shut off QCD processes because our description fails
  real(kind=rk), parameter, public   :: QCDcut = 0.6_rk
  real(kind=rk), parameter, public   :: omegax = 0.12_rk ! observed DM relic density
  real(kind=rk), parameter, public   :: rhocrit = 10.537_rk ! critical density, in GeV/m^3
  real(kind=rk), parameter, public   :: s0 = 2.8912e9_rk ! entropy today, in m^(-3)
  real(kind=rk), parameter, public   :: v = 246.0_rk ! Higgs vev in GeV
  real(kind=rk), parameter, public   :: mh = 125.0_rk ! Higgs mass in GeV


  ! number of rhs
  integer(kind=ik), public :: nrhs != 3

  contains
  !initialize global communicator of MPI_World
  subroutine set_mpi_comm_global(comm)
     implicit none
       integer, intent(in) :: comm

     BOLTZMANN_COMM=comm
 end subroutine set_mpi_comm_global

 subroutine abort_it(msg)
  implicit none
  character(len=*), intent(in) :: msg
  integer(kind=ik) :: mpierr

  write(*,*) msg
  !call abort
  call MPI_ABORT( BOLTZMANN_COMM, 666, mpierr)
end subroutine abort_it

end module
