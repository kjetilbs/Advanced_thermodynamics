####################################################################################################
# Plotting results
#
# Plotting results from the phase equilibrium calculation, using the constructed 
# ansArray from equilibriumCalculations.jl to extract the necessary data. The 
# following data is extracted, the plotted using pgfplots:
#   - (V,S) isotherm
#   - (S,H) isotherm
#   - (V,p) isotherm
#   - (T,S) isochor
#   - Liquid phase volume fraction
#   - Vapor phase volume fraction
#
# Author:     Kjetil Sonerud
# Updated:    2014-11-30 11:53:32
####################################################################################################

tic();

# Plotting on/off
plotting = false

##################################################
# Extracting isotherm data
##################################################
fh_VS_isoT = open("VS_isotherm.dat", "w")
fh_SH_isoT = open("SH_isotherm.dat", "w")
fh_VP_isoT = open("VP_isotherm.dat", "w")

row, col = size(ansArray[7])

for i = 1:3:row
    for j = 1:col
        # Writing (V,S)-data
        write(fh_VS_isoT, string(ansArray[3][i,j])*" "*string(ansArray[8][i,j])*"\n")
        # Writing (V,p)-data
        write(fh_VP_isoT, string(ansArray[3][i,j])*" "*string(ansArray[2][i,j])*"\n")
        # Writing (S,H)-data
        write(fh_SH_isoT, string(ansArray[8][i,j])*" "*string(ansArray[9][i,j])* "\n")
    end
    # Line break between data series
    write(fh_VS_isoT, "\n\n")
    write(fh_VP_isoT, "\n\n")
    write(fh_SH_isoT, "\n\n")
end

# Close files
close(fh_VS_isoT)
close(fh_VP_isoT)
close(fh_SH_isoT)

##################################################
# Extracting isochor data
##################################################
fh_TS_isoV = open("TS_isochor.dat", "w")

for i = 1:3:col
    for j = 1:row
        # Writing (T,S)-data
        write(fh_TS_isoV, string(ansArray[1][j,i])*" "*string(ansArray[8][j,i])*"\n")
    end
    # Line break between data series
    write(fh_TS_isoV, "\n\n")
end

# Close files
close(fh_TS_isoV)

##################################################
# Extracting volume fraction data
##################################################
# Opening files
fh_liq_phase    = open("liq_phase_composition.dat", "w")
fh_vap_phase    = open("vap_phase_composition.dat", "w")
fh_iterations   = open("iterations.dat", "w")

# Writing headers
write(fh_liq_phase, "V_tot T fraction \n")
write(fh_vap_phase, "V_tot T fraction \n")
write(fh_iterations, "V_tot T iterations \n")

# Dimensions of results
row, col = size(ansArray[7])

# Looping over the results
for i = 1:row
    for j = 1:col
        # Results stored in ansArray:
        #   - [1]: ansTemperature
        #   - [3]: ansVolumeTotal
        #   - [4]: ansVolumeVap
        #   - [5]: ansVolumeLiq

        # Writing liquid phase data
        write(fh_liq_phase, string(ansArray[3][i,j])*" "*string(ansArray[1][i,j])*" "*string(ansArray[5][i,j]/ansArray[3][i,j])*"\n")
        
        # Writing vapor phase data
        write(fh_vap_phase, string(ansArray[3][i,j])*" "*string(ansArray[1][i,j])*" "*string(ansArray[4][i,j]/ansArray[3][i,j])*"\n")
        
        # Writing number of iterations
        write(fh_iterations, string(ansArray[3][i,j])*" "*string(ansArray[1][i,j])*" "*string(ansArray[10][i,j])*"\n")
    end
    # Line break between data series
    # write(fh_liq_phase, "\n\n")
    # write(fh_vap_phase, "\n\n")
end

# Close files
close(fh_liq_phase)
close(fh_vap_phase)
close(fh_iterations)

if plotting == true
    ########################################
    # Plotting: (V,S)-diagram
    ########################################
    run(`pdflatex VS_isotherm.tex`)
    run(`pdflatex VS_isotherm.tex`)
    
    sleep(1)
    # Remove clutter
    run(`rm VS_isotherm.aux`)
    run(`rm VS_isotherm.log`)
    
    run(`open VS_isotherm.pdf`)
    
    ########################################
    # Plotting: (V,p)-diagram
    ########################################
    run(`pdflatex VP_isotherm.tex`)
    run(`pdflatex VP_isotherm.tex`)
    
    sleep(1)
    # Remove clutter
    run(`rm VP_isotherm.aux`)
    run(`rm VP_isotherm.log`)
    
    run(`open VP_isotherm.pdf`)
    
    ########################################
    # Plotting: (T,S)-diagram
    ########################################
    run(`pdflatex TS_isochor.tex`)
    run(`pdflatex TS_isochor.tex`)
    
    sleep(1)
    # Remove clutter
    run(`rm TS_isochor.aux`)
    run(`rm TS_isochor.log`)
    
    run(`open TS_isochor.pdf`)
    
    ########################################
    # Plotting: (S,H)-diagram
    ########################################
    run(`pdflatex SH_isotherm.tex`)
    run(`pdflatex SH_isotherm.tex`)
    
    sleep(1)
    # Remove clutter
    run(`rm SH_isotherm.aux`)
    run(`rm SH_isotherm.log`)
    
    run(`open SH_isotherm.pdf`)
    
    ########################################
    # Plotting: liquid phase diagram
    ########################################
    run(`pdflatex liq_phase_composition.tex`)
    
    sleep(1)
    # Remove clutter
    run(`rm liq_phase_composition.aux`)
    run(`rm liq_phase_composition.log`)
    
    run(`open liq_phase_composition.pdf`)
    
    ########################################
    # Plotting: vapor phase diagram
    ########################################
    run(`pdflatex vap_phase_composition.tex`)
    
    sleep(1)
    # Remove clutter
    run(`rm vap_phase_composition.aux`)
    run(`rm vap_phase_composition.log`)
    
    run(`open vap_phase_composition.pdf`)

    ########################################
    # Plotting: number of iterations
    ########################################
    run(`pdflatex iterations.tex`)
    
    sleep(1)
    # Remove clutter
    run(`rm iterations.aux`)
    run(`rm iterations.log`)
    
    run(`open iterations.pdf`)
end

toc();
