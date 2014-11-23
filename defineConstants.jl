#############################################################
# Define constants and read component data
#
# Define constants needed to perform the equilibrium calculations
# using the Redlich-Kwong equation of state. 
#
# Read component data from component_data.csv, which contains 
# component data (Tc,Pc,href,sref,cp(@175K)) for components
# N2, CH4, C2H6, C3H8, C4H10. The data is used for calculating 
# Redlich-Kwong parameters and thermodynamic properties of the
# component mixture. 
# 
# The data is downloaded from DIPPR:
# http://dippr.byu.edu/students/chemsearch.asp
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-20 13:04:33
#############################################################

# Read in the component data
componentData = readdlm("component_data.csv",',','\n',skipstart=2)

#############################################################
# Extracting data
#############################################################

# Component names
const global componentNames = convert(Array{String,1},componentData[:,1])

# Critical temperatures [kK]
const global tc = convert(Array{Float64,1},componentData[:,2])*1e-3

# Critical pressure [MPa]
const global pc = convert(Array{Float64,1},componentData[:,3])*1e-6

# Reference enthalpies [kJ/mol]
const global h_ref = convert(Array{Float64,1},componentData[:,4])*1e-6

# Reference entropies [kJ/kK mol]
const global s_ref = convert(Array{Float64,1},componentData[:,5])*1e-3

# Ideal gas heat capacity @175K [kJ/kK mol]
# Note: assumed to be constant. Should be changed for accuracy
const global c_p = convert(Array{Float64,1},componentData[:,6])*1e-3

# Cp-values; ideally, these should be integrated analytically 
# from the DIPPR expressions, and then used as is. As a quick
# fix, constant values may be assumed - they will not affect 
# equilibrium calculations themselves anyway

#############################################################
# Thermodynamic constants
#############################################################

# Universal gas constant [kJ/kK mol]
const global R = 8.3145*1e-3

# Reference temperature [kK]
const global T_ref = 298.15*1e-3

# Reference pressure [MPa]
const global p_ref = 1e5*1e-6

# Redlich-Kwong constant a [kJ^2 kK^0.5/kmol^2 kPa]
const global a_RK = (1/(9*((2^(1/3))-1)))*((R^2*tc.^(5/2))./pc)

# Redlich-Kwong constant b [kJ/kPa kmol]
const global b_RK = (((2^(1/3))-1)/3)*((R*tc)./pc)

#############################################################
# Other variables 
#############################################################

# Debug, level 1?
const global DEBUG1 = false

# Debug, level 2?
const global DEBUG2 = true