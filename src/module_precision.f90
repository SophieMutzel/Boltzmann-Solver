module module_precision
  use mpi
  implicit none
  integer, parameter, public         :: rk=selected_real_kind(8)
  integer, parameter, public         :: ik=selected_int_kind(8)

  real(kind=rk), parameter, public   :: pi  = 4.0 * atan(1.0)
  integer(kind=ik)                   :: BOLTZMANN_COMM

  contains
  !initialize global communicator of MPI_World
  subroutine set_mpi_comm_global(comm)
     implicit none
       integer, intent(in) :: comm

     BOLTZMANN_COMM=comm
 end subroutine set_mpi_comm_global
end module
