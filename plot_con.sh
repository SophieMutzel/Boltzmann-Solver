#!/bin/bash
gnuplot plot_rhs.gp
pdflatex plots/rhs.tex
mv rhs.* plots/
open -a Preview plots/rhs.pdf
