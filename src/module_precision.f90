module module_precision
  implicit none
  integer, parameter, public      :: rk=selected_real_kind(8)
  integer, parameter, public      :: ik=selected_int_kind(8)

  real(kind=rk), parameter, public   :: pi  = 4.0 * atan(1.0)
end module
