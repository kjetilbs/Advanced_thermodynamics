#!/bin/bash
# My bash script for doing Redlich-Kwong-plotting
# Written by: 	Kjetil Sonerud
# Updated: 		2014-11-16 17:42:23

# Run script with:  
#	./sinus_plot_test.sh 
# after using: 
#	chmod u+x sinus_plot_test.sh 
# to set the correct permissions

echo "Making figure..."

# Generating data using Julia
julia main.jl
# epstopdf introduction-inc.eps
#pdflatex -interaction=batchmode introduction.tex
cd plotting
pdflatex -halt-on-error pV_plot.tex | grep -a3 ^!
open pV_plot.pdf
cd ..

# Output
echo "Figure completed!"