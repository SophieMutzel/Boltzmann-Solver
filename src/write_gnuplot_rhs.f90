! write preface for gnuplot file which plots evolution of rhs terms
subroutine write_gnuplot_rhs(unit_number, filename)
  implicit none
  character(len=*), intent(in)      :: filename
  integer(kind=ik), intent(in)      :: unit_number
  integer(kind=ik)                  :: io_error

  open (unit=unit_number, file=filename//".gp", status='replace', action='write', iostat=io_error)

  if ( io_error == 0) then
    write(unit_number,*)  "set term tikz standalone color size 5in,3in"
    write(unit_number,*)  "set border linewidth 1.5"
    write(unit_number,*)  "set output '", filename,".tex'"
    write(unit_number,*)  "set xlabel '$\log_{10}(z)=\log_{10}(m_\chi/T)$'"
    write(unit_number,*)  "set ylabel 'rhs'"
    write(unit_number,*)  "set key font ',1'"
    write(unit_number,*)  "set logscale y"
    write(unit_number,*)  "set yrange[1e-20:1e-6]"
    write(unit_number,*)  "set grid ytics lc rgb '#bbbbbb' lw 1 lt 0"
    write(unit_number,*)  "set grid xtics lc rgb '#bbbbbb' lw 1 lt 0"
    write(unit_number,*)  "load 'gnuplot-palettes-master/dark2.pal'"
    write(unit_number,*)  "plot '", filename, ".txt' using 1:2 with lines ls   1 title '$H(T)$' , \"!'$T^\prime$' , \"
    write(unit_number,*)  "'", filename, ".txt' using 1:3 with lines ls   2 title '$\langle \sigma v_{\chi \bar\chi \leftrightarrow aa}\rangle n_\chi$', \"!'$\chi \chi \leftrightarrow aa$', \"
    write(unit_number,*)  "'", filename, ".txt' using 1:4 with lines ls   3 title '$\langle \sigma v_{aa \leftrightarrow \chi \bar\chi}\rangle n_a$', \"!'SM $\rightarrow \chi$', \"
    write(unit_number,*)  "'", filename, ".txt' using 1:5 with lines ls   4 title '$\langle \sigma v_{a\chi \leftrightarrow a \chi}\rangle n_a$', \"!'$3 H n_\chi$', \"
    write(unit_number,*)  "'", filename, ".txt' using 1:6 with lines ls   5 title '$\langle \sigma v_{\chi\bar\chi \leftrightarrow f \bar f}\rangle n_\chi$', \"!'$3 H n_a$', \"
    write(unit_number,*)  "'", filename, ".txt' using 1:6 with lines ls   5 dashtype 2 title '$\langle \sigma v_{\chi\bar\chi \leftrightarrow f \bar f}\rangle n_\chi$', \"!'$3 H n_a$', \"
    write(unit_number,*)  "'", filename, ".txt' using 1:7 with lines ls   6 title '$\langle \sigma v_{a SM \leftrightarrow SM SM}\rangle n_{SM}$', \"!'SM $\rightarrow a$', \"
    write(unit_number,*)  "'", filename, ".txt' using 1:7 with lines ls   6 dashtype 2 title '$\langle \sigma v_{a SM \leftrightarrow SM SM}\rangle n_{SM}$', \"!'SM $\rightarrow a$', \"
  else
    write(*,*) 'error', io_error,' while opening the file ', filename
  end if
  close(unit_number)
end subroutine write_gnuplot_rhs
