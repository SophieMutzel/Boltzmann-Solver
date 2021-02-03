set term tikz standalone color size 5in,3in
set border linewidth 1.5
set output 'plots/boltzmann.tex'
set xlabel '$\log_{10}(z)=\log_{10}(m_\chi/T)$'
set ylabel '$Y_\chi$'
set key font ',1'
set key outside
set logscale y
set yrange[1e-15:1e3]
set grid ytics lc rgb '#bbbbbb' lw 1 lt 0
set grid xtics lc rgb '#bbbbbb' lw 1 lt 0
load 'gnuplot-palettes-master/set1.pal'
plot 'temp/try.txt' using 1:2 smooth bezier ls   1 title '$Y_\chi(T)$' , \
      'temp/try.txt' using 1:3 smooth bezier ls   2 title '$Y_a(T)$' , \
      'temp/try.txt' using 1:4 smooth bezier ls   3 title '$Y_{\chi,eq}(T)$' , \
      'temp/try.txt' using 1:5 smooth bezier ls   4 title '$Y_{\chi,eq}(T^\prime)$' , \
      'temp/try6.txt' using 1:2 smooth bezier ls   7 title '$Y_\chi(T)$ with eq without con' , \
      'temp/try6.txt' using 1:3 smooth bezier ls   8 title '$Y_a(T)$ with eq without con', \
      'temp/try6.txt' using 1:6 smooth bezier ls   9 title '$Y_{a,eq}(T^\prime)$' , \
      'temp/try6.txt' using 1:7 smooth bezier ls   10 title '$z^\prime$' , \
      'temp/try7.txt' using 1:2 smooth bezier ls   11 title '$Y_\chi(T)$' , \
      'temp/try7.txt' using 1:3 smooth bezier ls   12 title '$Y_a(T)$' , \
