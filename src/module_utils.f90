module module_utils

  use module_precision
  use mpi

  implicit none
  type, public :: type_argsint
    real(kind=rk)                     :: T, mf, nc, mx, ma
  end type type_argsint

  contains
    subroutine linspace(from, to, array)

        implicit none

        real(kind=rk), intent(in)                  :: from, to
        real(kind=rk), dimension(:), intent(inout) :: array
        real(kind=rk)                              :: range
        integer(kind=ik)                           :: n, i

        n = size(array,1)
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from
            return
        end if

        do i=1, n
            array(i) = from + range * real(i - 1) / real(n - 1)
        end do

    end subroutine
! --------------------------------------------------------------------------!
! Took these routines from micromegas
! --------------------------------------------------------------------------!
    real(kind=rk) function bessI0( x )
      implicit none
      real(kind=rk), intent(in)   :: x
      real(kind=rk)               :: ax, ans
      real(kind=rk)               :: y

      ax = abs(x)
      if (ax < 3.75) then
        y = x/3.75
        y = y*y
        ans = (1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492 &
        +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))))
      else
        y=3.75/ax
        ans= ((exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1 &
        +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2 &
        +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1 &
        +y*0.392377e-2)))))))))
      end if
      bessI0 = ans
      return
    end function bessI0

    real(kind=rk) function bessI1( x )
      implicit none
      real(kind=rk), intent(in) :: x
      real(kind=rk)             :: ax, ans
      real(kind=rk)             :: y

      ax = abs(x)
      if (ax < 3.75) then
        y = x/3.75
        y = y*y
        ans=(ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934 &
        +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))))))
      else
        y=3.75/ax
        ans=(0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1 &
        -y*0.420059e-2)))
        ans= (0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2 &
        +y*(0.163801e-2+y*(-0.1031555e-1+y*ans)))))
        ans = ans * ((exp(ax)/sqrt(ax)))
      end if
      if (x < 0.0) then
        bessI1 = -ans
      else
        bessI1 = ans
      end if
      return
    end function bessI1

    real(kind=rk) function bessK0( x )
      !M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
      !Applied Mathematics Series vol. 55 (1964), Washington.
      implicit none
      real(kind=rk), intent(in) :: x
      real(kind=rk)             :: y, ans

      if (x <= 2.0) then
        y=x*x/4.0
        ans=(-log(x/2.0)*bessI0(x))+(-0.57721566+y*(0.42278420 &
        +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2 &
        +y*(0.10750e-3+y*0.74e-5))))))
      else
        y=2.0/x;
        ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1 &
        +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2 &
        +y*(-0.251540e-2+y*0.53208e-3))))))
      end if
      bessK0 = ans
      return
    end function bessK0

    real(kind=rk) function bessK1( x )
      implicit none
      real(kind=rk), intent(in) :: x
      real(kind=rk)             :: y, ans

      if (x <= 2.0) then
        y=x*x/4.0
        ans=(log(x/2.0)*bessI1(x))+(1.0/x)*(1.0+y*(0.15443144 &
        +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1 &
        +y*(-0.110404e-2+y*(-0.4686e-4)))))))
      else
        y=2.0/x
        ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619 &
        +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2 &
        +y*(0.325614e-2+y*(-0.68245e-3)))))))
      end if
      bessK1 = ans
      return
    end function bessK1

    real(kind=rk) function bessK2( x )
      implicit none
      real(kind=rk), intent(in) :: x
      real(kind=rk)             :: bk, bkm, bkp, tox

      tox = 2.0/x
      bkm = bessK0(x)
      bk  = bessK1(x)
      bkp = bkm+tox*bk
      bkm = bk
      bk  = bkp
      bessK2 = bk
      return
    end function bessK2

    integer(kind=ik) function str2int(str)
      implicit none
      character(len=*), intent(in)  :: str

      read(str,*)  str2int
      return
    end function str2int
    character(len=2) function int2str(int)
      implicit none
      integer(kind=ik), intent(in)  :: int

      write(int2str,'(I2)') int
      return
    end function int2str
    character(len=3) function float2str(float)
      implicit none
      real(kind=rk), intent(in)  :: float

      write(float2str,'(F1.1)') float
      return
    end function float2str

    include "quadpack.f90"
    include "interpolation.f90"
end module
