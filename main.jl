#############################################################
# Main.jl
#
# Main.jl is the master file for calculating phase equilibria
# for a multicomponent system using a minimum Helmholtz energy
# formulation to calculate the equlibrium. 
#  
# The components are: 
# N2, CH4, C2H6, C3H8, C4H10
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-06 21:24:02
#############################################################

# Clear all variables
workspace()

# Reading component data
include("readComponentData.jl")

# Defining constants
include("defineConstants.jl")

# Initialization of the parameters
n_feed = [0.016 0.9450 0.026 0.0081 0.0052];
n_feed_tot = sum(n_feed);

println("Sum feed: " * string(n_feed_tot));

# Important; find good starting values for gas and liquid
# when it comes to phase distribution

# The same is true for volume; initial values could be
# 	- Ideal gas volume for vapor phase
# 	- 1.1 * hard-sphere volume for liquid phase

# Constants
# R = 8.3145	# Universal gas constant, [J/K mol]

# # Critical parameters for the components
# # Species: N2, CH4, C2H6, C3H8, C4H10
# Tc 	= [126.2 190.564 ]; # [K]
# Pc 	= [3.4e6 4.599e6 ]; # [Pa]

# # Thermodynamic data for the components
# s_0 = [1.91609e5 1.86270e5] # [J/K mol]
# cp = [29.106 33.499]; #[J/K mol]

# Constant Cp-values calculated from DIPPR @175K: 

# Cp-values; ideally, these should be integrated analytically 
# from the DIPPR expressions, and then used as is. As a quick
# fix, constant values may be assumed - they will not affect 
# equilibrium calculations themselves anyway

