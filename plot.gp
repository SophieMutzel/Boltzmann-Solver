set term tikz standalone color size 5in,3in
set border linewidth 1.5
set output 'plots/trycon.tex'
set xlabel '$\log_{10}(z)=\log_{10}(m_\chi/T)$'
set ylabel '$Y_\chi$'
set key font ',1'
set logscale y
set yrange[1e-15:1e3]
set grid ytics lc rgb '#bbbbbb' lw 1 lt 0
set grid xtics lc rgb '#bbbbbb' lw 1 lt 0
load 'gnuplot-palettes-master/set1.pal'
plot 'temp/mx25.96E-11.63E-01.txt' using 1:2  with lines ls   1 dashtype 2 title '$Y_\chi(T)$' , \
      'temp/mx25.96E-11.63E-01.txt' using 1:3  with lines ls   2 dashtype 2 title '$Y_a(T)$' , \
      'temp/mx25.96E-11.63E-01.txt' using 1:4  with lines ls   3 title '$Y_{\chi,eq}(T)$' , \
      'temp/mx25.96E-11.63E-01.txt' using 1:5  with lines ls   1  title '$Y_{\chi,eq}(T^\prime)$' , \
      'temp/mx25.96E-11.63E-01.txt' using 1:6  with lines ls   2  title '$Y_{a,eq}(T^\prime)$' , \
      'temp/mx25.96E-11.63E-01.txt' using 1:7  with lines ls   6 title '$z^\prime$' , \
      8.18e-12   ls 7 title 'observed relic density'
#plot 'temp/compAoife.10E-10.10E-01.txt' using 1:2  with lines ls   1 dashtype 2 title '$Y_\chi(T)$' , \
#      'temp/compAoife.10E-10.10E-01.txt' using 1:3  with lines ls   2 dashtype 2 title '$Y_a(T)$' , \
#      'temp/compAoife.10E-10.10E-01.txt' using 1:4  with lines ls   3 title '$Y_{\chi,eq}(T)$' , \
#      'temp/compAoife.10E-10.10E-01.txt' using 1:5  with lines ls   1  title '$Y_{\chi,eq}(T^\prime)$' , \
#      'temp/compAoife.10E-10.10E-01.txt' using 1:6  with lines ls   2  title '$Y_{a,eq}(T^\prime)$' , \
#      'temp/compAoife.10E-10.10E-01.txt' using 1:7  with lines ls   6 title '$z^\prime$' , \
#      'temp/compAoife2.10E-10.10E-01.txt' using 1:2  with lines ls   1 dashtype 3 title '$Y_\chi(T)$' , \
#      'temp/compAoife2.10E-10.10E-01.txt' using 1:3  with lines ls   2 dashtype 3 title '$Y_a(T)$' , \
#      'temp/compAoifecon.10E-10.10E-01.txt' using 1:2  with lines ls   1 dashtype 4 title '$Y_\chi(T)$' , \
#      'temp/compAoifecon.10E-10.10E-01.txt' using 1:3  with lines ls   2 dashtype 4 title '$Y_a(T)$'
