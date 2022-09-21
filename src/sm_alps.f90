! a-SM interaction processes, kernel functions for integration of reaction rates
real(kind=rk) function kernel_afgf( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(in)     :: argsint
    kernel_afgf = sigma_afgf(s,argsint%mf,argsint%ma,argsint%g)*4.0_rk*&
                  F(s,argsint%ma,argsint%mf)*F(s,argsint%ma,argsint%mf)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
    return
end function kernel_afgf

real(kind=rk) function kernel_afgf_th( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(inout)  :: argsint
    kernel_afgf_th = sigma_afgf_th(s,argsint%mf,argsint%ma,argsint%mg,argsint%g,argsint)*4.0_rk*&
                  F(s,argsint%ma,argsint%mf)*F(s,argsint%ma,argsint%mf)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
    return
end function kernel_afgf_th

real(kind=rk) function kernel_agff( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(in)     :: argsint
    kernel_agff = sigma_agff(s,argsint%mf,argsint%nc,argsint%ma,argsint%g)*4.0_rk*&
                  F(s,argsint%ma,0.0_rk)*F(s,argsint%ma,0.0_rk)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
    return
end function kernel_agff

real(kind=rk) function kernel_agff_th( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(in)     :: argsint
    kernel_agff_th = sigma_agff_th(s,argsint%mf,argsint%nc,argsint%ma,argsint%mg,argsint%g)*4.0_rk*&
                  F(s,argsint%ma,argsint%mg)*F(s,argsint%ma,argsint%mg)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
    return
end function kernel_agff_th

real(kind=rk) function kernel_ahff( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(in)     :: argsint
    kernel_ahff = sigma_ahff(s,argsint%mf,argsint%nc,argsint%ma,argsint%g)*4.0_rk*&
                  F(s,argsint%ma,mh)*F(s,argsint%ma,mh)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
    return
end function kernel_ahff

real(kind=rk) function kernel_ffag( s, argsint )
  implicit none
    real(kind=rk), intent(in)           :: s
    type (type_argsint), intent(in)     :: argsint
    kernel_ffag = sigma_ffag(s,argsint%mf,argsint%nc,argsint%ma,argsint%g)*4.0_rk*&
                  F(s,argsint%mf,argsint%mf)*F(s,argsint%mf,argsint%mf)/sqrt(s) * bessK1( sqrt(s)/argsint%T)
    return
end function kernel_ffag
