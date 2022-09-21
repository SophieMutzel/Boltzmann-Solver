#!/bin/bash
mxx=10
maa=1
k=-13
x=-1.97511
filename=mx${mxx}ma${maa}e${k}e${x}
filerhs=rhs_mx${mxx}ma${maa}e${k}e${x}
./general --kappa ${k},${k} --gaxx ${x},${x}  --zstart -1 --zmax 2 --dzplot 0.02 --file $filename --mx ${mxx} --ma ${maa}
# --grid
gnuplot temp/${filename}.gp
pdflatex temp/${filename}.tex
rm ${filename}.aux
rm ${filename}.tex
rm ${filename}.log
mv ${filename}.pdf plots/
open -a Preview plots/${filename}.pdf
gnuplot temp/${filerhs}.gp
pdflatex temp/${filerhs}.tex
rm ${filerhs}.aux
rm ${filerhs}.tex
rm ${filerhs}.log
mv ${filerhs}.pdf plots/
open -a Preview plots/${filerhs}.pdf
