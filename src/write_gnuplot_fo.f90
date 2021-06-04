subroutine write_gnuplot_fo(unit_number, filename, Yobs)
  implicit none
  character(len=*), intent(in)      :: filename
  integer(kind=ik), intent(in)      :: unit_number
  real(kind=rk), intent(in)         :: Yobs
  integer(kind=ik)                  :: io_error

  open (unit=unit_number, file=filename//".gp", status='replace', action='write', iostat=io_error)

  if ( io_error == 0) then
    write(unit_number,*) "set term tikz standalone color size 5in,3in"
    write(unit_number,*) "set border linewidth 1.5"
    write(unit_number,*) "set output '", filename, ".tex'"
    write(unit_number,*) "set xlabel '$\log_{10}(z)=\log_{10}(m_\chi/T)$'"
    write(unit_number,*) "set ylabel '$Y_\chi$'"
    write(unit_number,*) "set key font ',1'"
    write(unit_number,*) "set logscale y"
    write(unit_number,*) "set yrange[1e-15:0.1]"
    write(unit_number,*) "set grid ytics lc rgb '#bbbbbb' lw 1 lt 0"
    write(unit_number,*) "set grid xtics lc rgb '#bbbbbb' lw 1 lt 0"
    write(unit_number,*) "load 'gnuplot-palettes-master/paired.pal'"
    write(unit_number,*) "plot '", filename,".txt' using 1:4  with lines ls   1 title '$Y_{\chi,eq}(T)$' , \"
    write(unit_number,*) "'", filename, ".txt' using 1:2  with lines ls   2 title '$Y_\chi(T)$' , \"
    write(unit_number,*) "'", filename,".txt' using 1:5  with lines ls  3  title '$Y_{a,eq}(T)$' , \"
    write(unit_number,*) "'", filename, ".txt' using 1:3  with lines ls   4 title '$Y_a(T)$' , \"
    write(unit_number,*) Yobs, "   ls 6 title 'observed relic density'"
  else
    write(*,*) 'error', io_error,' while opening the file ', filename
  end if
  close(unit_number)
end subroutine write_gnuplot_fo
