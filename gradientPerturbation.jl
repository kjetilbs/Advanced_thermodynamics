#############################################################
# Gradient perturbation
#
# Gradient perturbation (in an attempt) to check if the 
# Hessian is correct
#
# Author:     Kjetil Sonerud
# Updated:    2014-11-23 16:32:58
#############################################################

workspace()

#############################################################
# Import modules
#############################################################
# Defining constants
include("defineConstants.jl")

# Using ideal gas module
include("idealGas.jl")
using idealGas

# Using Redlich-Kwong module
include("redlichKwong.jl")
using redlichKwong

#############################################################
# Import modules
#############################################################
# Feed composition of natural gas, ref. Gonzalez & Lee (1968)
n_feed = [0.016, 0.9450, 0.026, 0.0081, 0.0052]

# Normalizing; converting to vector in Julia
n_feed = n_feed/sum(n_feed)

# Distribution of components in the two phases; light 
# components primarily in the gas 
n_vap = [0.99; 0.9; 0.5; 0.1; 0.01].*n_feed

n_liquid = n_feed - n_vap

# Temperature
T       = 0.180     # [kK]
V_tot   = 5e-4      # [m^3]
V_liq   = 1.5*redlichKwong.redlichKwongB(n_liquid)
V_vap   = V_tot - V_liq

#############################################################
# Calculate the gradient and Hessian
#############################################################
# Calculating the current gradient vector for the vapor phase
p_vap           = pressure(T,V_vap,n_vap)
mu_vap          = chemicalPotential(T,V_vap,n_vap)

# Calculating the current gradient vector for the vapor phase 
gradient_vap    = [p_vap, mu_vap]

# Calculating the current Hessian matrix for the vapor phase
H_vap           = hessian(T,V_vap,n_vap)

# Small volume perturbation
delta_V_vap     = 1e-10
V_vap           += delta_V_vap

# Re-calculating the gradient vector for the vapor phase 
p_vap           = pressure(T,V_vap,n_vap)
mu_vap          = chemicalPotential(T,V_vap,n_vap)

# Re-calculating gradient vector for the vapor phase 
perturbed_gradient_vap    = [p_vap, mu_vap]

# Storing the "Hessian" from perturbed gradient vectors
H_vap_from_grad = zeros(H_vap)

# Result of volume perturbation
H_vap_from_grad[1,:] = (gradient_vap - perturbed_gradient_vap)/delta_V_vap

# println(round(H_vap,4))
# println("\n")
# println(round(H_vap_from_grad,4))
# println("\n")

# Mole number perturbation
delta_n_vap = 1e-10

# "Reset" volume
V_vap           -= delta_V_vap

for i in 1:length(n_vap)
    # Perturbed mole vector
    perturbed_n_vap = n_vap

    # Small mole number perturbation
    perturbed_n_vap[i]       += delta_n_vap

    # Re-calculating the gradient vector for the vapor phase 
    p_vap           = pressure(T,V_vap,perturbed_n_vap)
    mu_vap          = chemicalPotential(T,V_vap,perturbed_n_vap)

    # Re-calculating gradient vector for the vapor phase 
    perturbed_gradient_vap    = [p_vap, mu_vap]

    # "Hessian" from perturbed gradient
    H_vap_from_grad[i+1,:] = (gradient_vap - perturbed_gradient_vap)/delta_n_vap
    # println(string(i)*" : "*string(round(H_vap_from_grad,4)))
end

#############################################################
# Printing results
#############################################################

println("Actual Hessian")
println(round(H_vap,4))
println("\n")
println("Numeric Hessian")
println(round(H_vap_from_grad,4))
println("\n")

# Calculating the relative error
println("Relative error comparing actual and numeric Hessian")
println(round(((H_vap-H_vap_from_grad)./H_vap),3))

