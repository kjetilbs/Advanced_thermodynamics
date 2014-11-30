################################################################################
# Define constants and read component data
#
# Define constants needed to perform the equilibrium calculations using the 
# Redlich-Kwong equation of state. 
#
# Read component data from component_data.csv, which contains component data 
# $T_c$, $P_c$, $h_{\mathrm{ref}}$, $s_{\mathrm{ref}}$, $c_p$ (@ 175$K$) 
# for components N2, CH4, C2H6, C3H8, C4H10. 
#
# The data is used for calculating Redlich-Kwong parameters and thermodynamic 
# properties of the component mixture. 
# 
# The data is downloaded from DIPPR:
# http://dippr.byu.edu/students/chemsearch.asp
#
# Author:   Kjetil Sonerud
# Updated:  2014-11-26 16:09:50
################################################################################

# Read in the component data
componentData = readdlm("component_data.csv",',','\n',skipstart=2)

################################################################################
# Extracting data
################################################################################

# Component names
const global componentNames = convert(Array{String,1},componentData[:,1])

# Critical temperatures $\mathrm{[K]}$
const global tc = convert(Array{Float64,1},componentData[:,2])

# Critical pressure $\mathrm{[Pa]}$
const global pc = convert(Array{Float64,1},componentData[:,3])

# Reference enthalpies $\mathrm{[J/mol]}$
const global h_ref = convert(Array{Float64,1},componentData[:,4])*1e-3

# Reference entropies $\mathrm{[J/K mol]}$
const global s_ref = convert(Array{Float64,1},componentData[:,5])*1e-3

# Ideal gas heat capacity @175$K$ $\mathrm{[J/K mol]}$
const global c_p = convert(Array{Float64,1},componentData[:,6])*1e-3

# Regarding the $c_p$-values: as these will not affect equilibrium calculations 
# they are assumed to be constant and independent of $T$. This is an approximation, 
# and ideally, the $c_p$-expressions as functions of temperature should be used.

################################################################################
# Thermodynamic constants
################################################################################

# Universal gas constant $\mathrm{[J/K mol]}$
const global R      = 8.3145

# Reference temperature $\mathrm{[K]}$
const global T_ref  = 298.15

# Reference pressure $\mathrm{[Pa]}$
const global p_ref  = 1e5

# Redlich-Kwong constant a $\mathrm{[kJ^2 kK^{0.5}/kmol^2 kPa]}$
const global a_RK   = (1/(9*((2^(1/3))-1)))*((R^2*tc.^(5/2))./pc)

# Redlich-Kwong constant b $\mathrm{[kJ/kPa kmol]}$
const global b_RK   = (((2^(1/3))-1)/3)*((R*tc)./pc)

################################################################################
# Other variables 
################################################################################

# Turn debugging on/off
const global DEBUG = false
