################################################################################
# (p,V)-plot using Redlich-Kwong
#
# (p,V)-plot using Redlich-Kwong equation of state. n is fixed, 
# T varied in outer loop, V varied in inner loop. The resulting 
# output is stored and plotted externally using pgfplots
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-30 17:54:15
################################################################################

################################################################################
# Phase composition
################################################################################

# Feed composition of natural gas, ref. Gonzalez & Lee (1968)
n_feed = [0.016 0.9450 0.026 0.0081 0.0052];

# Normalizing; converting to vector in Julia
n_feed = vec(n_feed/sum(n_feed))

# Distribution of components in the two phases
# n_vapor = [0.9; 0.8; 0.1; 0.01; 0.001].*n_feed
n_vap = [0.9; 0.7; 0.3; 0.01; 0.01].*n_feed

n_liq = n_feed - n_vap

################################################################################
# Redlich-Kwong plotting; $(V_i,p_i)$-diagram
################################################################################

include("redlichKwong.jl")
using redlichKwong

# Defining filenames
filename_vap = "sandbox/vap_Vp_data.dat"
filename_liq = "sandbox/liq_Vp_data.dat"
# Open files to write results to	
f_handle_vap = open(filename_vap, "w")
f_handle_liq = open(filename_liq, "w")

# Total data
filename = "results/only_feed_Vp_data.dat"
f_handle = open(filename, "w")

# Outer loop: temperature
for T in 0.150:0.005:0.230
	# Inner loop; volume
	for V in logspace(log10(1.1*redlichKwong.redlichKwongB(n_liq)),0,100)
		# Calculate $p_{vap}$ from RK EOS
		p_vap = pressure(T,V,n_vap)
		write(f_handle_vap, string(V) * " " * string(1e6*p_vap) * "\n")
		
		# # Calculate $p_{liq}$ from RK EOS
		p_liq = pressure(T,V,n_liq)
		write(f_handle_liq, string(V) * " " * string(1e6*p_liq) * "\n")

		# Calculate $p$ from RK EOS
		p = pressure(T,V,n_feed)
		write(f_handle, string(V) * " " * string(1e6*p) * "\n")
	end
	# Separate data series
	write(f_handle_vap, "\n\n")
	write(f_handle_liq, "\n\n")
	write(f_handle, "\n\n")
end

# Close files
close(f_handle_vap)
close(f_handle_liq)
close(f_handle)

################################################################################
# Plotting
################################################################################
# cd("plotting")
# run(`pdflatex tot_pV_plot.tex`)

# # Remove clutter
# run(`rm tot_pV_plot.aux`)
# run(`rm tot_pV_plot.log`)

# cd("..")
# run(`open plotting/tot_pV_plot.pdf`)



