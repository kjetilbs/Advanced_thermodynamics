################################################################################
# (V,p)-plot using Redlich-Kwong, plotting particular solutions
#
# (V,p)-plot using Redlich-Kwong equation of state. n is fixed, T varied in outer
# loop, V varied in inner loop. The resulting output is stored and plotted 
# externally using pgfplots. The point of the plot is to verify that the actual
# solutions lie on the (V,p) curve for a particular isotherm
#
# Author:   Kjetil Sonerud
# Updated:  2014-11-16 15:02:33
################################################################################

################################################################################
# Phase composition
################################################################################

# Feed composition of natural gas, ref. Gonzalez & Lee (1968)
n_feed = [0.016 0.9450 0.026 0.0081 0.0052];

# Normalizing; converting to vector in Julia
n_feed = vec(n_feed/sum(n_feed))

################################################################################
# Redlich-Kwong plotting; $(V_i,p_i)$-diagram
################################################################################

include("redlichKwong.jl")
using redlichKwong

# Outer loop: temperature
for T in [170,180,190]
    # File name 
    filename = "plotting/Vp_data_for_feed_composition_"*string(T)*"K.dat"
    f_handle = open(filename, "w")

    # Inner loop; volume
    for V in logspace(log10(1.1*redlichKwong.redlichKwongB(n_feed)),0,100)
        # Calculate $p$ from RK EOS
        p = pressure(T,V,n_feed)
        write(f_handle, string(V) * " " * string(p) * "\n")
    end
    # Close files
    close(f_handle)
end

################################################################################
# Plotting
################################################################################
cd("plotting")
run(`pdflatex tot_Vp_isotherms_with_results.tex`)

# Remove clutter
run(`rm tot_Vp_isotherms_with_results.aux`)
run(`rm tot_Vp_isotherms_with_results.log`)

cd("..")
run(`open plotting/tot_Vp_isotherms_with_results.pdf`)



