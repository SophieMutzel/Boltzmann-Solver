set term tikz standalone color size 5in,3in
set border linewidth 1.5
set output 'plots/con.tex'
set xlabel '$\log_{10}(z)=\log_{10}(m_\chi/T)$'
set ylabel 'rhs'
set key font ',1'
set key outside
set logscale y
set yrange[1e-20:1e-8]
set grid ytics lc rgb '#bbbbbb' lw 1 lt 0
set grid xtics lc rgb '#bbbbbb' lw 1 lt 0
load 'gnuplot-palettes-master/dark2.pal'
plot 'temp/con.txt' using 1:2 smooth bezier ls   1 title '$\chi$:  $aa\rightarrow \chi \chi$' , \
      'temp/con.txt' using 1:3 smooth bezier ls   2 title '$\chi$:  $\chi \chi \rightarrow aa$' , \
      'temp/con.txt' using 1:4 smooth bezier ls   3 title '$\chi$:  $f f  \rightarrow \chi \chi$' , \
      'temp/con.txt' using 1:5 smooth bezier ls   4 title '$ a$:  $ff \rightarrow a i$' , \
      'temp/con2.txt' using 1:2 smooth bezier ls  5 title '$\chi$:  $aa\rightarrow \chi \chi$' , \
      'temp/con2.txt' using 1:3 smooth bezier ls  6 title '$\chi$:  $\chi \chi \rightarrow a a$', \
      'temp/con3.txt' using 1:2 smooth bezier ls  7 title '$\chi$:  $aa\rightarrow \chi \chi$ eq' , \
      'temp/con3.txt' using 1:3 smooth bezier ls  8 title '$\chi$:  $\chi \chi \rightarrow a a$ eq', \
      'temp/con3.txt' using 1:4 smooth bezier ls  9 title '$H(z)$'
