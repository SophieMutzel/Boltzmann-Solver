#!/bin/bash
mxx=10
maa=1
#for k in -15 -14 -13 -12 -11
#do
#k=( -15 -14 -13 -12 -11 )
#x=( -2.38782 -2.16869 -1.95032 -1.72374 -1.4995 )
#k=( -12 )
#x=( -1.72374 )
k=( -16 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1)
#k=( -10 -9 -8 )
x=( -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 )
#  for x in -1.6 -1.8 -2 -2.2 -2.4
#  do
for (( n=0; n<=14; n++ ))
  do
    for (( i=0; i<=10; i++ ))
      do
        filename=mx${mxx}ma${maa}e${k[$n]}e${x[$i]}
        filerhs=rhs_mx${mxx}ma${maa}e${k[$n]}e${x[$i]}
    #./boltzmann --kappa ${k},${k} --gaxx ${x},${x}  --zstart -1.4 --zmax 3 --dzplot 0.05 --file mx10e${k}e${x} --mx 10
        ./general --kappa ${k[$n]},${k[$n]} --gaxx ${x[$i]},${x[$i]}  --zstart -1 --zmax 0 --dzplot 0.05 --file $filename --mx ${mxx} --ma ${maa}
#    gnuplot temp/${filename}.gp
#    pdflatex temp/${filename}.tex
#    rm ${filename}.aux
#    rm ${filename}.tex
#    rm ${filename}.log
#    mv ${filename}.pdf plots/
#    open -a Preview plots/${filename}.pdf
#    gnuplot temp/${filerhs}.gp
#    pdflatex temp/${filerhs}.tex
#    rm ${filerhs}.aux
#    rm ${filerhs}.tex
#    rm ${filerhs}.log
#    mv ${filerhs}.pdf plots/
#    open -a Preview plots/${filerhs}.pdf
  done
done
