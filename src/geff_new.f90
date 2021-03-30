real(kind=rk) function geff_rho_h(T)
  implicit none
  real(kind=rk), intent(in)     :: T
  real(kind=rk), dimension(12)  :: a = (/1.0_rk,1.117240_rk,3.12672E-01_rk,&
                  -4.68049E-02_rk, -2.65004E-02_rk, -1.19760E-03_rk,&
                  1.82812E-04_rk,  1.36436E-04_rk,  8.55051E-05_rk,&
                  1.22840E-05_rk,  3.82259E-07_rk,  -6.87035E-09_rk/)

  real(kind=rk), dimension(12)  :: b = (/1.43382E-02_rk,1.37559E-02_rk,2.92108E-03_rk,&
                  -5.38533E-04_rk, -1.62496E-04_rk, -2.87906E-05_rk,&
                  -3.84278E-06_rk, 2.78776E-06_rk,  7.40342E-07_rk,&
                  1.17210E-07_rk,  3.72499E-09_rk,  -6.74107E-11_rk/)
  real(kind=rk)                 :: LT, x, temp1, temp2
  integer(kind=ik)              :: i

  LT = log(T)
  x = 1.0_rk
  temp1 = 0.0_rk
  temp2 = 0.0_rk

  do i = 1, 12
    temp1 = temp1 + a(i) * x
    temp2 = temp2 + b(i) * x
    x = x * LT
  end do

  geff_rho_h = temp1 / temp2
  return
end function geff_rho_h

real(kind=rk) function ratio(T)
  implicit none
  real(kind=rk), intent(in)     :: T
  real(kind=rk), dimension(12)  :: a = (/1.0_rk,6.07869E-01_rk,-1.54485E-01_rk,&
                  -2.24034E-01_rk, -2.82147E-02_rk, 2.90620E-02_rk,&
                  6.86778E-03_rk,  -1.05E-03_rk, -1.69104E-04_rk,&
                  1.06301E-05_rk,  1.69528E-06_rk,  -9.33311E-08_rk/)
  real(kind=rk), dimension(12)  :: b = (/7.07388E+01_rk,9.18011E+01_rk,3.31892E+01_rk,&
                  -1.39779E+00_rk, -1.52558E+00_rk, -1.97857E-02_rk,&
                  -1.60146E-01_rk, 8.22615E-05_rk,  2.02651E-02_rk,&
                  -1.82134E-05_rk, 7.83943E-05_rk,  7.13518E-05_rk/)
  real(kind=rk)                 :: LT, x, temp1, temp2
  integer(kind=ik)              :: i
  LT = log(T)
  x = 1.0_rk
  temp1 = 0.0_rk
  temp2 = 0.0_rk
  do i=1,12
    temp1 = temp1 + a(i) * x
    temp2 = temp2 + b(i) * x
    x = x * LT
  end do
  ratio = temp1 / temp2 + 1.0_rk
  return
end function ratio

real(kind=rk) function geff_s_h(T)
  implicit none
  real(kind=rk), intent(in)     :: T
  geff_s_h = 1.0_rk / ratio(T) * geff_rho_h(T)
  return
end function geff_s_h
real(kind=rk) function fr(x)
  implicit none
  real(kind=rk), intent(in)     :: x
  fr =  exp(-1.04855_rk * x) * (1.0_rk + 1.03757_rk * x + &
        0.508630_rk * x * x + 0.0893988_rk * x * x * x)
  return
end function fr

real(kind=rk) function br(x)
  implicit none
  real(kind=rk), intent(in)     :: x
  br = exp(-1.03149_rk * x) * (1.0_rk + 1.03317_rk * x +&
        0.398264_rk * x * x + 0.0648056_rk * x * x * x)
  return
end function br

real(kind=rk) function fs(x)
  implicit none
  real(kind=rk), intent(in)     :: x
  fs =  exp(-1.04190_rk * x) * (1.0_rk + 1.03400_rk * x &
        + 0.456426_rk * x * x + 0.0595248_rk * x * x * x)
  return
end function fs

real(kind=rk) function bs(x)
  implicit none
  real(kind=rk), intent(in)     :: x
  bs = exp(-1.03365_rk * x) * (1.0_rk + 1.03397_rk * x + &
      0.342548_rk * x * x + 0.0506182_rk * x * x * x)
  return
end function bs

real(kind=rk) function Sfit(x)
  implicit none
  real(kind=rk), intent(in)     :: x
  Sfit = 1.0_rk + 7.0_rk / 4.0_rk * exp(-1.0419_rk * x) * &
        (1.0_rk + 1.034_rk * x + 0.456426_rk * x * x + 0.0595249_rk * x * x * x)
  return
end function Sfit

real(kind=rk) function geff_rho_l(T) !=1,=1
  implicit none
  real(kind=rk), intent(in)   :: T
  real(kind=rk)               :: me = 511e-6_rk, mmu = 0.1056_rk, mpi0 = 0.135_rk
  real(kind=rk)               :: mpip = 0.140_rk, m1 = 0.5_rk, m2 = 0.77_rk, m3 = 1.2_rk, m4 = 2.0_rk

  geff_rho_l = 2.030_rk + 1.353_rk * Sfit(me/T)**(4.0_rk/3.0_rk) + 3.495_rk * fr(me / T) +&
         3.446_rk * fr(mmu / T) + 1.05_rk * br(mpi0 / T) + 2.08_rk * br(mpip / T) +&
         4.165_rk * br(m1 / T) + 30.55_rk * br(m2 / T) + 89.4_rk * br(m3 / T) +&
         8209_rk * br(m4 / T)
  return
end function geff_rho_l

real(kind=rk) function geff_s_l(T)
  implicit none
  real(kind=rk), intent(in)   :: T
  real(kind=rk)               :: me = 511e-6_rk, mmu = 0.1056_rk, mpi0 = 0.135_rk
  real(kind=rk)               :: mpip = 0.140_rk, m1 = 0.5_rk, m2 = 0.77_rk, m3 = 1.2_rk, m4 = 2.0_rk

  geff_s_l = 2.008_rk + 1.923_rk * Sfit(me / T) + 3.442_rk * fs(me / T) +&
         3.468_rk * fs(mmu / T) + 1.034_rk * bs(mpi0 / T) + 2.068_rk * bs(mpip / T) +&
         4.160_rk * bs(m1 / T) + 30.55_rk * bs(m2 / T) + 90.0_rk * bs( m3 / T) +&
         6209_rk * br(m4 / T) ! p1*181.7*bs(p2*m3/T)
  return
end function geff_s_l

real(kind=rk) function geff_rho(T)
  implicit none
  real(kind=rk), intent(in)     :: T
  if (T < 0.12_rk) then
    geff_rho = geff_rho_l(T)
  else
    geff_rho = geff_rho_h(T)
  end if
  return
end function geff_rho
real(kind=rk) function geff_s(T)
  implicit none
  real(kind=rk), intent(in)     :: T
  if (T < 0.12_rk) then
    geff_s = geff_s_l(T)
  else
    geff_s = geff_s_h(T)
  end if
  return
end function geff_s
