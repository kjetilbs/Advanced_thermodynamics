#############################################################
# (p,V)-plot using Redlich-Kwong
#
# (p,V)-plot using Redlich-Kwong equation of state. n is fixed, 
# T varied in outer loop, V varied in inner loop. The resulting 
# output is stored and plotted externally using pgfplots
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-16 15:02:33
#############################################################

# Outer loop: temperature
for T in 175:5:210
	# Open file to write results to
	filename = "plotting/pVdata_T_"*string(T)*".dat"
	f_handle = open(filename, "w")
	# Testing 
	# write(f_handle, "Partyhatt! \n")
	println(redlichKwongParameterA(T,n_vapor))
	for V in logspace(log10(1.1*redlichKwongParameterB(n_vapor)),0,100)
		# Calculate p from RK EOS
		p = redlichKwongEOS(T,V,n_vapor)
		# println(p)
		write(f_handle, string(V) * " " * string(p) * "\n")
	end
	# Close file
	close(f_handle)
end

