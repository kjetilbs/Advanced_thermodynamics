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

# # Reading component data
# include("readComponentData.jl")

# # Defining constants
# include("defineConstants.jl")

# # Initialization of the parameters
# n_feed = [1 0 0 0 0]

# # To obtain correct dimensions 
# n_feed = vec(n_feed)

# Feed composition of natural gas, ref. Gonzalez & Lee (1968)
n_feed = [0.016 0.9450 0.026 0.0081 0.0052];

# Normalizing; converting to vector in Julia
n_feed = vec(n_feed/sum(n_feed))

# Distribution of components in the two phases
# n_vapor = [0.9; 0.8; 0.1; 0.01; 0.001].*n_feed
n_vapor = [0.5; 0.5; 0.5; 0.5; 0.5].*n_feed

n_liquid = n_feed - n_vapor

# println(n_vapor)
# println(n_liquid)

# n_feed_tot = sum(n_feed);

# println("Sum feed: " * string(n_feed_tot));

# Important; find good starting values for gas and liquid
# when it comes to phase distribution

# The same is true for volume; initial values could be
# 	- Ideal gas volume for vapor phase
# 	- 1.1 * hard-sphere volume for liquid phase

# println(n_feed)

# println(n_vapor)
# println(n_vapor*2)
# println(n_feed)

# # Redlich-Kwong functions
# include("redlichKwong.jl")

# # T = 298.15
# # V = 0.1

# # p = redlichKwongEOS(T,V,n_feed)

# include("redlichKwongPlotting.jl")
# include("idealGas.jl")

# Importing ideal gas module
include("idealGasFunctions.jl")
importall idealGas
println("Contents of ideal gas module:")
whos(idealGas)
println("\n")

# Importing Redlich-Kwong module
include("redlichKwong.jl")
importall redlichKwong
println("Contents of Redlich-Kwong module:")
whos(redlichKwong)

# println(idealHessian(298.15,0.1,n_vapor))
# println(idealGasEOS(298.15,0.1,n_vapor))




