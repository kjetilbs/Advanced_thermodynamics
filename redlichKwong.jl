#############################################################
# Redlich-Kwong functions
#
# Functions related to the Redlich-Kwong equation of state
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-16 14:10:33
#############################################################

#############################################################
# Thermodynamic functions using Redlich-Kwong
#
# Functions for calculating thermodynamic potentials, 
# properties and useful derivatives using Redlich-Kwong EOS: 
#
#	- idealGasH: enthalpy for an ideal gas
#	- idealGasS: entropy for an ideal gas
#	- idealGasMu: chemical potential for an ideal gas
#	- idealGasA: Helmholtz free energy for an ideal gas
#
# Derivatives of Helmholtz free energy using ideal gas:
#	- Aig_T: 	$\pdc{A\ig}{T}{V,\vt{n}} = - S\ig$ 
#	- Aig_V: 	$\pdc{A\ig}{V}{T,\vt{n}} = - p\ig$
#	- Aig_n: 	$\pdc{A\ig}{\vt{n}}{T,V} = \mu\ig$
#	- Aig_TT:	$\pddc{A\ig}{T}{T}{V,\vt{n}}$ 
#	- Aig_TV:	$\pddc{A\ig}{T}{V}{\vt{n}}$ 
#	- Aig_Tn:	$\pddc{A\ig}{T}{\vt{n}}{V}$
#	- Aig_VV:	$\pddc{A\ig}{V}{V}{T,\vt{n}}$
#	- Aig_nV:	$\pddc{A\ig}{\vt{n}}{V}{T}$
#	- Aig_nn:	$\pddc{A\ig}{\vt{n}}{\vt{n}}{T,V}$
#
# Resulting Hessian matrix based on Helmholtz free energy
# for ideal gas with constant temperature $T$: 
# 	- idealGasHessian 
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-20 11:28:33
#############################################################

module redlichKwong

# Functions to be available in global scope once the module is
# imported/being used
export redlichKwongEOS							#rkHessian

# Reading component data
include("readComponentData.jl")

# Defining constants
include("defineConstants.jl")

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

#############################################################
# End module
#############################################################
end