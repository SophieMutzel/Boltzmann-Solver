#!/bin/bash
gnuplot plot2.gp
pdflatex plots/boltzmann2.tex
mv boltzmann2.* plots/
open -a Preview plots/boltzmann2.pdf
