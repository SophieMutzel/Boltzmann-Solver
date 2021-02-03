#!/bin/bash
gnuplot plot_con.gp
pdflatex plots/con.tex
mv con.* plots/
open -a Preview plots/con.pdf
