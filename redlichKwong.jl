#############################################################
# Redlich-Kwong functions
#
# Functions related to the Redlich-Kwong equation of state
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-16 14:10:33
#############################################################

#############################################################
# Redlich-Kwong EOS
#############################################################
function redlichKwongParameterA(T,n)
	# Calculate the parameter A for the Redlich-Kwong EOS
	# 
	#	A(T,n) = (n'*(a*a')*n)/(T^(0.5))
	# 
	# where a is defined in defineConstants.jl 

	# Using "dot" for dot product in Julia
	A = dot(n,(((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)))
end

function redlichKwongParameterB(n)
	# Calculate the parameter B for the Redlich-Kwong EOS
	# 
	# 	B(n) = b'*n
	#
	# 	where b is defined in defineConstants.jl  

	# Using "dot" for dot product in Julia
	B = dot(b_RK,n)
end

function redlichKwongEOS(T,V,n)
	# Calculate the pressure from the Redlich-Kwong EOS:
	# 	p_RK(T,V,n) = (NRT)/(V-b) - (A)/(V(V-B))
	
	# Calculating A and B parameters
	A = redlichKwongParameterA(T,n)
	B = redlichKwongParameterB(n)

	# Calculating the pressure
	p = (sum(n)*R*T)/(V-B) - (A/(V*(V+B)))
end