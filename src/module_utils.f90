! useful other routines, like Bessel functions etc.
module module_utils

  use module_precision
  use mpi

  implicit none
  type, public :: type_argsint
    real(kind=rk)                     :: T, mf, nc, mx, ma, g, mg, ra_ini, s
    real(kind=rk), allocatable        :: drhoa(:,:)
    logical                           :: helper
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
    !Returns the modified Bessel function Kn(x) for positive x and n ≥ 2.
    real (kind=rk) FUNCTION bessk_s(n,x)
      IMPLICIT NONE
      INTEGER(kind=ik), INTENT(IN)  :: n
      REAL(kind=rk), INTENT(IN)     :: x
      INTEGER(kind=ik)              :: j
      REAL(kind=rk)                 :: bk,bkm,bkp,tox
      !call assert(n >= 2, x > 0.0, ’bessk_s args’)
      tox=2.0_rk/x
      bkm=bessk0(x)
      bk=bessk1(x)
      do j=1,n-1
        bkp=bkm+j*tox*bk
        bkm=bk
        bk=bkp
      end do
      bessk_s=bk
    END FUNCTION bessk_s

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

    character(len=7) function exp2str(float)
      implicit none
      real(kind=rk), intent(in)  :: float

      write(exp2str,'(E7.2)') float
      return
    end function exp2str
    function inv(A) result(Ainv)
      real(kind=rk), dimension(:,:), intent(in) :: A
      real(kind=rk), dimension(size(A,1),size(A,2)) :: Ainv

      real(kind=rk), dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info

      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
    end function inv
    subroutine get_bezier_coeff(points, n, A, B)
      implicit none
      real(kind=rk), intent(in)     :: points(n+1,2)
      integer(kind=ik), intent(in)  :: n
      real(kind=rk), intent(out)    :: A(n,2), B(n,2)
      real(kind=rk), dimension(n,n) :: C
      real(kind=rk), dimension(n,2) :: P
      integer(kind=ik)              :: i

      !n = len(points) - 1
      C=0.0_rk
      C(1,1) = 2.0_rk
      do i=2,n
        C(i,i) = 4.0_rk
        C(i-1,i) = 1.0_rk
        C(i,i-1) = 1.0_rk
      end do
      C(n,n) = 7.0_rk
      C(n,n-1) = 2.0_rk
      do i = 1, n
        P(i,:) = 2.0_rk * ( 2.0_rk * points(i,:) + points(i+1,:) )
      end do
      P(1,:) = points(1,:) + 2.0_rk * points(2,:)
      P(n,:) = 8.0_rk * points(n,:) + points(n+1,:)

      A(:,1) = matmul(inv(C),P(:,1))
      A(:,2) = matmul(inv(C),P(:,2))
      do i = 1, n-1
        B(i,:) = 2.0_rk * points(i+1,:) - A(i+1,:)
      end do
      B(n,:) = (A(n,:) + points(n+1,:))/2.0_rk

    end subroutine get_bezier_coeff

    subroutine get_cubic(a,b,c,d,t,curve)
      !returns the general Bezier cubic formula given 4 control points
      implicit none
      real(kind=rk), intent(in), dimension(2)   :: a,b,c,d
      real(kind=rk), intent(in)                 :: t
      real(kind=rk), intent(out),dimension(2)   :: curve(2)
      curve(:) = (1.0_rk-t)**3 * a + 3.0_rk * (1.0_rk-t)**2 * t * b &
                  + 3.0_rk*(1.0_rk-t)*t*t*c + t*t*t*d
    end subroutine get_cubic

    ! return one cubic curve for each consecutive points
!    subroutine get_bezier_cubic(points, n, bezier)
!      implicit none
!      real(kind=rk), intent(in)     :: points(n+1,2)
!      integer(kind=ik), intent(in)  :: n
!      real(kind=rk), intent(out)    :: bezier(n,2)
!      real(kind=rk)                 :: A(n,2), B(n,2)
!      integer(kind=ik)              :: i
!
!      call get_bezier_coeff(points,n,A,B)
!      do i=1,n
!        call get_cubic(points(i,:),A(i,:),B(i,:),points(i+1,:),t,bezier(i,:))
!      end do
!
!    end subroutine get_bezier_cubic
    subroutine test_bezier(x,y,x_eval,y_eval,n,A,B)
      implicit none
      real(kind=rk), intent(in)     :: x(:),y(:),x_eval,A(:,:),B(:,:)
      real(kind=rk), intent(out)    :: y_eval
      integer(kind=ik), intent(in)  :: n
      real(kind=rk)                 :: xmin,ymin,xmax,ymax
      integer(kind=4)              :: l, j, i
      real(kind=rk)                 :: tvals
      real(kind=rk)                 :: points(n+1,2), bezier(2)!,A(n,2), B(n,2)

      xmin = minval(x)
      xmax = maxval(x)
      ymin = minval(y)
      ymax = maxval(y)
      points(:,1)=(x-xmin)/(xmax-xmin)
      points(:,2)=(y-ymin)/(ymax-ymin)
      !call linspace(0.0_rk,1.0_rk,tvals)
      !call get_bezier_coeff(points,n,A,B)
      call r8vec_bracket ( n+1, x, x_eval, i, j )
      tvals = (x_eval-x(i))/(x(i+1)-x(i))
      !l=1
          call get_cubic(points(i,:),A(i,:),B(i,:),points(i+1,:),tvals,bezier(:))
          !l = l+1
        !end do
      !end do
      !x_eval = bezier(1)*(xmax-xmin)+xmin
      y_eval = bezier(2)*(ymax-ymin)+ymin

    end subroutine test_bezier

    real(kind=rk) function alpha_s(mu)
      implicit none
      real(kind=rk), intent(in)     :: mu
      real(kind=rk)                 :: nf, Lambda, lmulambda, beta0, beta1, beta2

      if (mu<1.67_rk) then
        nf = 3.0_rk
        Lambda = 0.372_rk
      else
        if (mu<4.8_rk) then
          nf = 4.0_rk
          Lambda = 0.325_rk
        else
          if (mu<173_rk) then
            nf = 5.0_rk
            Lambda = 0.226_rk
          else
            nf = 6.0_rk
            Lambda = 0.092_rk
          end if
        end if
      end if
      lmulambda = log(mu*mu/Lambda/Lambda)
      beta0 = 11.0_rk-2.0_rk*nf/3.0_rk
      beta1 = 102.0_rk-38.0_rk*nf/3.0_rk
      beta2 = 2857.0_rk/2.0_rk-5033.0_rk*nf/18.0_rk+325.0_rk*nf*nf/54.0_rk
      alpha_s = 4.0_rk*pi/(beta0*lmulambda)*&
                (1.0_rk-beta1/beta0/beta0*log(lmulambda)/&
                lmulambda+beta1*beta1/(beta0*beta0*beta0*beta0*&
                lmulambda*lmulambda)*((log(lmulambda)-0.5_rk)*(log(lmulambda)-0.5_rk)&
                +beta2*beta0/beta1/beta1-5.0_rk/4.0_rk))
    end function alpha_s

    real(kind=rk) function alpha_qed_th(mu)
      implicit none
      real(kind=rk), intent(in)     :: mu
      real(kind=rk)                 :: alphaMZ, MZ, beta0

      beta0=-4.0_rk/3.0_rk
      MZ = 91.1876_rk
      alphaMZ = 1.0_rk/128.962_rk
      alpha_qed_th = alphaMZ/(1.0_rk+2.0_rk*beta0*alphaMZ/(4.0_rk*pi)*Log(mu/MZ))

    end function alpha_qed_th
!    subroutine test_bezier(x,y,x_eval,y_eval,n,nvals)
!      implicit none
!      real(kind=rk), intent(in)     :: x(:),y(:)
!      real(kind=rk), intent(out)    :: x_eval(:),y_eval(:)
!      integer(kind=ik), intent(in)  :: n,nvals
!      real(kind=rk)                 :: xmin,ymin,xmax,ymax
!      integer(kind=rk)              :: l, j, i
!      real(kind=rk)                 :: tvals(nvals)
!      real(kind=rk)                 :: A(n,2), B(n,2), points(n+1,2), bezier(n*nvals,2)
!
!      xmin = minval(x)
!      xmax = maxval(x)
!      ymin = minval(y)
!      ymax = maxval(y)
!      points(:,1)=(x-xmin)/(xmax-xmin)
!      points(:,2)=(y-ymin)/(ymax-ymin)
!      call linspace(0.0_rk,1.0_rk,tvals)
!      call get_bezier_coeff(points,n,A,B)
!      call r8vec_bracket ( n, x, t, left, right )
!      l=1
!      do i=1,n
!        do j=1,nvals
!          call get_cubic(points(i,:),A(i,:),B(i,:),points(i+1,:),tvals(j),bezier(l,:))
!          l = l+1
!        end do
!      end do
!      x_eval(:) = bezier(:,1)*(xmax-xmin)+xmin
!      y_eval(:) = bezier(:,2)*(ymax-ymin)+ymin
!
!    end subroutine test_bezier

    include "quadpack.f90"
    include "spline.f90"
    include "interpolation.f90"

end module
