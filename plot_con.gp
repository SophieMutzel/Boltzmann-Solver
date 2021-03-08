set term tikz standalone color size 5in,3in
set border linewidth 1.5
set output 'plots/gammas4.tex'
set xlabel '$\log_{10}(z)=\log_{10}(m_\chi/T)$'
set ylabel 'rhs'
set key font ',1'
set key outside
set logscale y
set yrange[1e-20:1e0]
set grid ytics lc rgb '#bbbbbb' lw 1 lt 0
set grid xtics lc rgb '#bbbbbb' lw 1 lt 0
load 'gnuplot-palettes-master/dark2.pal'
plot 'temp/con.txt' using 1:2 with lines ls   1 title '$\chi \chi \rightarrow aa$' , \
      'temp/con.txt' using 1:3 with lines ls   2 title '$a a \rightarrow \chi \chi$' , \
      'temp/con.txt' using 1:4 with lines ls   3 title '$ a$:  $ff \rightarrow a i$' , \
      'temp/con.txt' using 1:5 with lines ls   4 title '$ a$:  $ff \rightarrow \chi \chi$' , \
      'temp/con.txt' using 1:6 with lines ls  5 title '$H(z)$'
