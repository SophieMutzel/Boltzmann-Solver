subroutine ikv ( v, x, bk, bi )

!*****************************************************************************80
!
!! IKV compute modified Bessel function Iv(x) and Kv(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    17 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V, the order of Iv(x) and Kv(x).
!    V = N + V0.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) VM, the highest order computed.
!
!    Output, real ( kind = 8 ) BI(0:N), DI(0:N), BK(0:N), DK(0:N), the
!    values of In+v0(x), In+v0'(x), Kn+v0(x), Kn+v0'(x).
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) bi(1:*)
  real ( kind = 8 ) bi0
  real ( kind = 8 ) bk(1:*)
  real ( kind = 8 ) bk0
  real ( kind = 8 ) bk1
  real ( kind = 8 ) bk2
  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  real ( kind = 8 ) cs
  real ( kind = 8 ) ct
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) gan
  real ( kind = 8 ) gap
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  integer ( kind = 4 ) m
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) n
  real ( kind = 8 ) pi
  real ( kind = 8 ) piv
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) sum
  real ( kind = 8 ) v
  real ( kind = 8 ) v0
  real ( kind = 8 ) v0n
  real ( kind = 8 ) v0p
  real ( kind = 8 ) vm
  real ( kind = 8 ) vt
  real ( kind = 8 ) w0
  real ( kind = 8 ) wa
  real ( kind = 8 ) ww
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  pi = 3.141592653589793D+00
  x2 = x * x
  n = int ( v )
  v0 = v - n
  if ( n == 0 ) then
    n = 1
  end if

  if ( x < 1.0D-100 ) then

    do k = 0, n
      bi(k) = 0.0D+00
      bk(k) = -1.0D+300
    end do

    if ( v == 0.0D+00 ) then
      bi(0) = 1.0D+00
    end if

    vm = v
    return

  end if

  piv = pi * v0
  vt = 4.0D+00 * v0 * v0

  if ( v0 == 0.0D+00 ) then
    a1 = 1.0D+00
  else
    v0p = 1.0D+00 + v0
    call gamma ( v0p, gap )
    a1 = ( 0.5D+00 * x ) ** v0 / gap
  end if

  if ( x < 35.0D+00 ) then
    k0 = 14
  else if ( x < 50.0D+00 ) then
    k0 = 10
  else
    k0 = 8
  end if

  if ( x <= 18.0D+00 ) then

    bi0 = 1.0D+00
    r = 1.0D+00
    do k = 1, 30
      r = 0.25D+00 * r * x2 / ( k * ( k + v0 ) )
      bi0 = bi0 + r
      if ( abs ( r / bi0 ) < 1.0D-15 ) then
        exit
      end if
    end do

    bi0 = bi0 * a1

  else

    ca = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
    sum = 1.0D+00
    r = 1.0D+00
    do k = 1, k0
      r = -0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
      sum = sum + r
    end do
    bi0 = ca * sum

  end if

  m = msta1 ( x, 200 )

  if ( m < n ) then
    n = m
  else
    m = msta2 ( x, n, 15 )
  end if

  f2 = 0.0D+00
  f1 = 1.0D-100
  do k = m, 0, -1
    f = 2.0D+00 * ( v0 + k + 1.0D+00 ) / x * f1 + f2
    if ( k <= n ) then
      bi(k) = f
    end if
    f2 = f1
    f1 = f
  end do

  cs = bi0 / f
  do k = 0, n
    bi(k) = cs * bi(k)
  end do

  if ( x <= 9.0D+00 ) then

    if ( v0 == 0.0D+00 ) then

      ct = - log ( 0.5D+00 * x ) - 0.5772156649015329D+00
      cs = 0.0D+00
      w0 = 0.0D+00
      r = 1.0D+00
      do k = 1, 50
        w0 = w0 + 1.0D+00 / k
        r = 0.25D+00 * r / ( k * k ) * x2
        cs = cs + r * ( w0 + ct )
        wa = abs ( cs )
        if ( abs ( ( wa - ww ) / wa ) < 1.0D-15 ) then
          exit
        end if
        ww = wa
      end do

      bk0 = ct + cs

    else

      v0n = 1.0D+00 - v0
      call gamma ( v0n, gan )
      a2 = 1.0D+00 / ( gan * ( 0.5D+00 * x ) ** v0 )
      a1 = ( 0.5D+00 * x ) ** v0 / gap
      sum = a2 - a1
      r1 = 1.0D+00
      r2 = 1.0D+00
      do k = 1, 120
        r1 = 0.25D+00 * r1 * x2 / ( k * ( k - v0 ) )
        r2 = 0.25D+00 * r2 * x2 / ( k * ( k + v0 ) )
        sum = sum + a2 * r1 - a1 * r2
        wa = abs ( sum )
        if ( abs ( ( wa - ww ) / wa ) < 1.0D-15 ) then
          exit
        end if
        ww = wa
      end do

      bk0 = 0.5D+00 * pi * sum / sin ( piv )

    end if

  else

    cb = exp ( - x ) * sqrt ( 0.5D+00 * pi / x )
    sum = 1.0D+00
    r = 1.0D+00
    do k = 1, k0
      r = 0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
      sum = sum + r
    end do
    bk0 = cb * sum

  end if

  bk1 = ( 1.0D+00 / x - bi(1) * bk0 ) / bi(0)
  bk(0) = bk0
  bk(1) = bk1
  do k = 2, n
    bk2 = 2.0D+00 * ( v0 + k - 1.0D+00 ) / x * bk1 + bk0
    bk(k) = bk2
    bk0 = bk1
    bk1 = bk2
  end do

  vm = n + v0

  return
end
