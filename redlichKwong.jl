#############################################################
# Thermodynamic functions using Redlich-Kwong
#
# Redlich-Kwong equation of state:
#	- redlichKwongA:	Parameter A used in the R-K EOS
#	- redlichKwongB:	Parameter B used in the R-K EOS
#	- redlichKwongP:	pressure from the R-K EOS
#
# Functions for calculating thermodynamic potentials, 
# properties and useful derivatives using Redlich-Kwong EOS: 
#
# Residual functions using R-K EOS: 
#	- residualA: 		Residual Helmholtz free energy
#	- residualS:		Residual entropy
#	- residualMu:		Residual chemical potential
#
# Derivatives of Helmholtz free energy using ideal gas:
#	- Aig_T: 			$\pdc{A\ig}{T}{V,\vt{n}} = - S\ig$ 
#	- Aig_V: 			$\pdc{A\ig}{V}{T,\vt{n}} = - p\ig$
#	- Aig_n: 			$\pdc{A\ig}{\vt{n}}{T,V} = \mu\ig$
#	- Aig_TT:			$\pddc{A\ig}{T}{T}{V,\vt{n}}$ 
#	- Aig_TV:			$\pddc{A\ig}{T}{V}{\vt{n}}$ 
#	- Aig_Tn:			$\pddc{A\ig}{T}{\vt{n}}{V}$
#	- Aig_VV:			$\pddc{A\ig}{V}{V}{T,\vt{n}}$
#	- Aig_nV:			$\pddc{A\ig}{\vt{n}}{V}{T}$
#	- Aig_nn:			$\pddc{A\ig}{\vt{n}}{\vt{n}}{T,V}$
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
export redlichKwongP, residualA							#rkHessian

# Defining constants and reading component data. Needs to be 
# included within the module scope
include("defineConstants.jl")

# Using the ideal gas module
# include("idealGas.jl")
using idealGas
println("Contents of ideal gas module:")
whos(idealGas)
println("\n")

#############################################################
# Redlich-Kwong EOS
#############################################################
function redlichKwongA(T,n)
	# Calculate the parameter A for the Redlich-Kwong EOS
	# $A(T,\vt{n}) = \frac{1}{T^{1/2}}\vtt{n}(\vt{a}\:\vtt{a})\vt{n}$
	# where a is defined in defineConstants.jl as 
	# $a_i = \frac{1}{9(2^{1/3}-1)}\frac{R^2T_{c,i}^{5/2}}{p_{c,i}} = \Omega_a \frac{R^2T_{c,i}^{5/2}}{p_{c,i}}$
	A = dot(n,(((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)))
end

function redlichKwongB(n)
	# Calculate the parameter B for the Redlich-Kwong EOS
	# $B(\vt{n}) = \vtt{b}\:\vt{n}$
	# where b is defined in defineConstants.jl as
	# $b_i = \frac{2^{1/3}-1}{3}\frac{RT_{c,i}}{p_{c,i}} = \Omega_b \frac{RT_{c,i}}{p_{c,i}}$ 
	B = dot(b_RK,n)
end

function redlichKwongP(T,V,n)
	# Calculate the pressure from the Redlich-Kwong EOS:
	# $p\rk = (T,V,\vt{n}) = \frac{NRT}{V-B} - \frac{A}{V(V+B)}$
	
	# Calculating A and B parameters
	A = redlichKwongA(T,n)
	B = redlichKwongB(n)

	# Calculating the pressure
	p = (sum(n)*R*T)/(V-B) - (A/(V*(V+B)))
end

#############################################################
# Residual thermodynamic potentials using Redlich-Kwong EOS
#############################################################
function residualA(T,V,n)
	# Calculate the residual Helmholtz free energy using the
	# R-K EOS
	# $A\rrk(T,V,\vt{n}) = NRT\ln \left(\frac{V}{V-B}\right)+\frac{A}{B}\ln \left(\frac{V}{V+B}\right)$
	# Calculating A and B parameters
	A = redlichKwongParameterA(T,n)
	B = redlichKwongParameterB(n)

	helm = sum(n)*R*T*log(V/(V-B)) + (A/B)*log(V/(V+B))
end

#############################################################
# End module
#############################################################
end