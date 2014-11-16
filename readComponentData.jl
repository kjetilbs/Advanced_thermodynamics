#############################################################
# Read component data
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
# Updated:	2014-11-16 13:22:46
#############################################################

# Read in the component data
componentData = readdlm("component_data.csv",',','\n',skipstart=2)

#############################################################
# Extracting data
#############################################################

# Component names
const global componentNames = convert(Array{String,1},componentData[:,1])

# Critical temperatures [K]
const global tc = convert(Array{Float64,1},componentData[:,2])

# Critical pressure [Pa]
const global pc = convert(Array{Float64,1},componentData[:,3])

# Reference enthalpies [J/mol]
const global h_ref = convert(Array{Float64,1},componentData[:,4])/1000

# Reference entropies [J/K mol]
const global s_ref = convert(Array{Float64,1},componentData[:,5])/1000

# Ideal gas heat capacity @175K [J/K mol]
# Note: assumed to be constant. Should be changed for accuracy
const global c_p = convert(Array{Float64,1},componentData[:,6])/1000

