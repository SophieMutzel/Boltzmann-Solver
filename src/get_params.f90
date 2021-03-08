subroutine get_params(params)
  use f90getopt
  implicit none
  type(type_params), intent(inout) :: params
  type(option_s)                :: opts(13)
  character(len=20)             :: gaxx, kappa
  logical                       :: kin, gin

  opts(1) = option_s( "mx", .true., 'x' )
  opts(2) = option_s( "ma",  .true.,  'a' )
  opts(3) = option_s( "zmax",  .true.,  'm' )
  opts(4) = option_s( "zstart",  .true.,  's' )
  opts(6) = option_s( "dzplot",  .true.,  'p' )
  opts(7) = option_s( "N",  .true.,  'n' )
  opts(8) = option_s( "kappa",  .true.,  'k' )
  opts(9) = option_s( "gaxx",  .true.,  'g' )
  opts(10) = option_s( "nt",  .true.,  't' )
  opts(11) = option_s( "dz",  .true.,  'z' )
  opts(12) = option_s( "file",  .true.,  'f' )
  opts(13) = option_s( "help", .false.,  'h')

  ! If no options were committed
  if (command_argument_count() .eq. 0 ) then
    write (*,*) "Available Options: --mx -x --ma -a --zmax -m --zstart -s&
                --dzplot -p --N -n --kappa -k --gaxx -g -h --help"
  end if

  ! Process options one by one
  do
    select case( getopt( "h:x:a:m:s:p:n:k:g:t:z:f", opts ) ) ! opts is optional (for longopts only)
      case( char(0) )
        exit
      case( 'x' )
        read(optarg,*) params%mx
      case( 'a' )
        read(optarg,*) params%ma
      case( 'm' )
        read(optarg,*) params%z_max
      case( 's' )
        read(optarg,*) params%z_start
      case( 'p' )
        read(optarg,*) params%dz_plot
      case( 'n' )
        read(optarg,*) params%N_tot
      case( 't' )
        read(optarg,*) params%nt
      case( 'z' )
        read(optarg,*) params%dz
      case( 'k' )
        kin = .true.
        kappa = optarg
      case( 'g' )
        gin = .true.
        gaxx = optarg
      case( 'f' )
        params%file = optarg
      case( 'h' )
        write (*,*) "Available Options:  --mx -x --ma -a --zmax -m --zstart -s&
                    --dzplot -p --N -n --kappa -k --gaxx -g -h --help"
    end select
  end do
!read(memstring(10:len_trim(memstring)-2),* )
  if ( kin) then
      read(kappa(1:INDEX(kappa, ",")-1),*) params%kappa_range(1)
      read(kappa(INDEX(kappa, ",")+1:), *) params%kappa_range(2)
  end if
  if ( gin) then
    read(gaxx(1:INDEX(gaxx, ",")-1),*) params%gaxx_range(1)
    read(gaxx(INDEX(gaxx, ",")+1:), *) params%gaxx_range(2)
  end if

end subroutine get_params
