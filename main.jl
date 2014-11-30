################################################################################
# Main.jl
#
# Main.jl is the master file for calculating phase equilibria
# for a multicomponent system using a minimum Helmholtz energy
# formulation to calculate the equlibrium. 
#  
# The components that constitute the system are: 
# N2, CH4, C2H6, C3H8, C4H10
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-30 12:19:10
################################################################################

# Clear all variables
workspace()

# Timing script
tic();

# Defining constants
include("defineConstants.jl")

# Using ideal gas module
include("idealGas.jl")
using idealGas

# Using Redlich-Kwong module
include("redlichKwong.jl")
using redlichKwong

##################################################
# Initial guess and iteration grid
##################################################

# Feed composition of natural gas, ref. Gonzalez & Lee (1968)
n_feed      = [0.016 0.9450 0.026 0.0081 0.0052];

# Normalizing; converting to vector in Julia
n_feed      = vec(n_feed/sum(n_feed))

# Guessing a distribution of components in the vapor phase
n_vapor     = [0.9; 0.8; 0.3; 0.05; 0.1].*n_feed

n_liquid    = n_feed - n_vapor

# Iteration range, temperature and volume
gridSize    = 50
maxT        = 200
minT        = 160
maxV        = 1e-3
minV        = 3e-4

rangeT = linspace(minT, maxT, gridSize)
rangeV = linspace(minV, maxV, gridSize)

# Calculating the hard sphere volume
V_liq = 1.3*redlichKwong.redlichKwongB(n_liquid)
V_vap = maxV - V_liq

# Initial guess vector
#   - x_guess[1]:       Vapor phase volume 
#   - x_guess[2:end]:   Vapor phase mole vector

x_guess = [V_vap, n_vapor]

##################################################
# Calculating phase equilibrium
##################################################
# Using equilibrium module
include("equilibriumCalculations.jl")
using equilibriumCalculations

ansArray = phaseEquilibrium(x_guess,n_feed,rangeT,rangeV);

##################################################
# Extract and plot results
##################################################
# cd("results/")
# pwd()
# println("Plotting...")
# sleep(2)
# include("results/results.jl")
# println("Plotting done!")
# cd("..")

toc();




