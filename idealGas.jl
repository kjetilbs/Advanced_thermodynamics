#############################################################
# Ideal gas law functions
#
# Functions related to the ideal gas law equation of state
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-19 17:01:39
#############################################################

function idealGasEOS(T,V,n)
	# Calculate the pressure from the ideal gas EOS:
	# 	p_ig(T,V,n) = NRT/V
	p = (sum(n)*R*T)/V
end


