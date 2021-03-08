#!/bin/bash
name=trycon
gnuplot plot.gp
pdflatex plots/${name}.tex
mv ${name}.* plots/
open -a Preview plots/${name}.pdf
