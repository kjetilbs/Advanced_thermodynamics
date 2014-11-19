#############################################################
# Define constants
#
# Define constants needed to perform the equilibrium calculations
# using the Redlich-Kwong equation of state
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-16 13:47:35
#############################################################

#############################################################
# Thermodynamic constants
#############################################################

# Universal gas constant [J/K mol]
const global R = 8.3145

# Reference temperature [K]
const global T_ref = 298.15

# Reference pressure [Pa]
const global p_ref = 1e5

# Redlich-Kwong constant a [J^2 K^0.5/mol^2 Pa]
const global a_RK = (1/(9*((2^(1/3))-1)))*((R^2*tc.^(5/2))./pc)

# Redlich-Kwong constant b [J/Pa mol]
const global b_RK = (((2^(1/3))-1)/3)*((R*tc)./pc)