################################################################################
# Gradient perturbation
#
# Gradient perturbation (in an attempt) to check if the 
# Hessian is correct
#
# Author:     Kjetil Sonerud
# Updated:    2014-11-23 16:32:58
################################################################################

workspace()

################################################################################
# Import modules
################################################################################
# Defining constants
include("defineConstants.jl")

# Using ideal gas module
include("idealGas.jl")
using idealGas

# Using Redlich-Kwong module
include("redlichKwong.jl")
using redlichKwong

################################################################################
# Define variables (T,V,n)
################################################################################
# Feed composition of natural gas, ref. Gonzalez & Lee (1968)
n_feed = [0.016, 0.9450, 0.026, 0.0081, 0.0052]

# Normalizing; converting to vector in Julia
n_feed = n_feed/sum(n_feed)

# Distribution of components in the two phases; light 
# components primarily in the gas 
# n_vap = [0.99; 0.9; 0.5; 0.1; 0.01].*n_feed
n_vap = [0.9; 0.8; 0.3; 0.05; 0.05].*n_feed

n_liquid = n_feed - n_vap

# Temperature
T       = 0.180     # [kK]
V_tot   = 5e-4      # [m^3]
V_liq   = 1.5*redlichKwong.redlichKwongB(n_liquid)
V_vap   = V_tot - V_liq

################################################################################
# Central difference scheme
################################################################################
function centralDifference(f_forward, f_backward, h)
    # Using central difference scheme to estimate numeric
    # derivative using
    # $f'(x) \approx \frac{(f(x+0.5h) - f(x-0.5h))}{h}$
    dfun = (f_forward - f_backward)/h
end

################################################################################
# Calculate the gradient and Hessian
################################################################################
# Calculating the current gradient vector for the vapor phase
p_vap               = pressure(T,V_vap,n_vap)
mu_vap              = chemicalPotential(T,V_vap,n_vap)

# Calculating the current gradient vector for the vapor phase 
gradient_vap        = [p_vap, mu_vap]

# Calculating the current ideal gas Hessian matrix for the 
# vapor phase
idealH_vap          = idealHessian(T,V_vap,n_vap)

# Calculating the current residual Hessian matrix for the 
# vapor phase
resH_vap            = redlichKwong.residualHessian(T,V_vap,n_vap)

# Calculating the current Hessian matrix for the vapor phase
H_vap               = hessian(T,V_vap,n_vap)

# Small volume perturbation
delta_V_vap         = 1e-8
# Small mole number perturbation
delta_n_vap         = 1e-8

# Preparing for central difference
forward_V_vap       = V_vap + 0.5*delta_V_vap
backward_V_vap      = V_vap - 0.5*delta_V_vap

# Initializing the numeric Hessian from perturbed gradient vectors
numericH_vap = zeros(H_vap)

forward_gradient_vap    =   [
                                -pressure(T,forward_V_vap,n_vap), 
                                chemicalPotential(T,forward_V_vap,n_vap)
                            ]
backward_gradient_vap    =  [
                                -pressure(T,backward_V_vap,n_vap), 
                                chemicalPotential(T,backward_V_vap,n_vap)
                            ]

# Result of volume perturbation
numericH_vap[1,:]       = centralDifference(forward_gradient_vap, backward_gradient_vap, delta_V_vap)

# println(round(H_vap,4))
# println("\n")
# println(round(numericH_vap,4))
# println("\n")


for i in 1:length(n_vap)
    # Preparing for central difference
    forward_n_vap       = deepcopy(n_vap)
    backward_n_vap      = deepcopy(n_vap)

    forward_n_vap[i]    += 0.5*delta_n_vap
    backward_n_vap[i]   -= 0.5*delta_n_vap
    
    forward_gradient_vap    =   [
                                    -pressure(T,V_vap,forward_n_vap), 
                                    chemicalPotential(T,V_vap,forward_n_vap)
                                ]
    backward_gradient_vap    =  [
                                    -pressure(T,V_vap,backward_n_vap), 
                                    chemicalPotential(T,V_vap,backward_n_vap)
                                ]

    # Checking mole vectors
    println("Before perturbation:   "*string(n_vap))
    println("Forward perturbation:  "*string(forward_n_vap))
    println("Backward perturbation: "*string(backward_n_vap))
    println("\n")

    # # Re-calculating the gradient vector for the vapor phase 
    # p_vap           = pressure(T,V_vap,perturbed_n_vap)
    # mu_vap          = chemicalPotential(T,V_vap,perturbed_n_vap)

    # # Re-calculating gradient vector for the vapor phase 
    # perturbed_gradient_vap    = [p_vap, mu_vap]

    # "Hessian" from perturbed gradient
    numericH_vap[i+1,:] = centralDifference(forward_gradient_vap, backward_gradient_vap, delta_n_vap)
    # println(string(i)*" : "*string(round(numericH_vap,4)))
end

################################################################################
# Printing results
################################################################################
println("Ideal Hessian")
println(round(idealH_vap,5))
println("Residual Hessian")
println(round(resH_vap,5))
println("Actual Hessian (for R-K)")
println(round(H_vap,5))
println("Numeric Hessian")
println(round(numericH_vap,5))
println("\n")

# println("a_i constants")
# println(round(a_RK,8))
# println("b_i constants")
# println(round(b_RK,8))
# println("\n")

# Calculating the relative error
println("Relative error comparing actual and numeric Hessian")
println(round(((numericH_vap-H_vap)./numericH_vap),10))
println("Norm:")
println(round(norm((numericH_vap-H_vap)./numericH_vap),10))

println("\n")
println("Analytic Hessian transposed:")
println(((H_vap-H_vap')./H_vap))
println("Norm:")
println(round(norm((H_vap-H_vap')./H_vap),6))

################################################################################
# Investigating Helmholtz free energy; mole number perturbation
################################################################################
helmholtzFunction   = redlichKwong.residualHelmholtz(T,V_vap,n_vap)
# helmholtzEuler      = (pressure(T,V_vap,n_vap)-idealGas.idealPressure(T,V_vap,n_vap))*V_vap + dot(redlichKwong.residualChemicalPotential(T,V_vap,n_vap),n_vap)

# println("\nDifference in residual Helmholtz:")
# println(helmholtzFunction-helmholtzEuler) 

# Initializing
numeric_Mu_vap = zeros(n_vap)

for i in 1:length(n_vap)
    # Preparing for central difference
    forward_n_vap       = deepcopy(n_vap)
    backward_n_vap      = deepcopy(n_vap)

    forward_n_vap[i]   += 0.5*delta_n_vap
    backward_n_vap[i]  -= 0.5*delta_n_vap
    
    forward_A_vap       = redlichKwong.residualHelmholtz(T,V_vap,forward_n_vap)
    backward_A_vap      = redlichKwong.residualHelmholtz(T,V_vap,backward_n_vap)
    
    # (pressure(T,V_vap,forward_n_vap) - idealGas.idealPressure(T,V_vap,forward_n_vap))*V_vap + dot(redlichKwong.residualChemicalPotential(T,V_vap,forward_n_vap),forward_n_vap)

    # Checking mole vectors
    println("Before perturbation:   "*string(n_vap))
    println("Forward perturbation:  "*string(forward_n_vap))
    println("Backward perturbation: "*string(backward_n_vap))
    println("\n")

    # Result of mole number perturbation
    numeric_Mu_vap[i]  = centralDifference(forward_A_vap, backward_A_vap, delta_n_vap)
end

println("\nDifference in chemical potential:")
println("Chemical potential from function:")
println(redlichKwong.residualChemicalPotential(T,V_vap,n_vap))
println("Ideal gas chemical potential:")
println(idealGas.idealChemicalPotential(T,V_vap,n_vap))
println("Chemical potential from Helmholtz perturbation:")
println(numeric_Mu_vap)
println("Relative error:")
println((redlichKwong.residualChemicalPotential(T,V_vap,n_vap) - numeric_Mu_vap)./redlichKwong.residualChemicalPotential(T,V_vap,n_vap))
println("\n\n")

################################################################################
# Investigating Helmholtz free energy; temperature perturbation
################################################################################
helmholtzFunction   = redlichKwong.residualHelmholtz(T,V_vap,n_vap)
# helmholtzEuler      = (pressure(T,V_vap,n_vap)-idealGas.idealPressure(T,V_vap,n_vap))*V_vap + dot(redlichKwong.residualChemicalPotential(T,V_vap,n_vap),n_vap)

# println("\nDifference in residual Helmholtz:")
# println(helmholtzFunction-helmholtzEuler)
# println("\n")

# Small temperature perturbation
delta_T_vap         = 1e-8

# Preparing for central difference
forward_T_vap       = deepcopy(T)
backward_T_vap      = deepcopy(T)

forward_T_vap      += 0.5*delta_T_vap
backward_T_vap     -= 0.5*delta_T_vap

forward_A_vap       = redlichKwong.residualHelmholtz(forward_T_vap,V_vap,n_vap)
# pressure(forward_T_vap,V_vap,n_vap)*V_vap + dot(chemicalPotential(forward_T_vap,V_vap,n_vap),n_vap)
backward_A_vap      = redlichKwong.residualHelmholtz(backward_T_vap,V_vap,n_vap)


# Checking temperature
println("Before perturbation (T):   "*string(T))
println("Forward perturbation (T):  "*string(forward_T_vap))
println("Backward perturbation (T): "*string(backward_T_vap))
println("\n")

# # Checking pressure
# println("Before perturbation (p):   "*string(pressure(T,V_vap,n_vap)-idealGas.idealPressure(T,V_vap,n_vap)))
# println("Forward perturbation (p):  "*string(pressure(forward_T_vap,V_vap,n_vap)-idealGas.idealPressure(forward_T_vap,V_vap,n_vap)))
# println("Backward perturbation (p): "*string(pressure(backward_T_vap,V_vap,n_vap)-idealGas.idealPressure(backward_T_vap,V_vap,n_vap)))
# println("\n")

# # Checking chemical potential
# println("Before perturbation (mu):   "*string(redlichKwong.residualChemicalPotential(T,V_vap,n_vap)))
# println("Forward perturbation (mu):  "*string(redlichKwong.residualChemicalPotential(forward_T_vap,V_vap,n_vap)))
# println("Backward perturbation (mu): "*string(redlichKwong.residualChemicalPotential(backward_T_vap,V_vap,n_vap)))
# println("\n")

# Result of mole number perturbation
numeric_S_vap   = -centralDifference(forward_A_vap, backward_A_vap, delta_T_vap)


println("\nDifference in entropy:")
println("Entropy from function:")
println(redlichKwong.residualEntropy(T,V_vap,n_vap))
println("Entropy from Helmholtz perturbation:")
println(numeric_S_vap)
println("Relative error:")
println((redlichKwong.residualEntropy(T,V_vap,n_vap) - numeric_S_vap)/redlichKwong.residualEntropy(T,V_vap,n_vap))

################################################################################
# Investigating Helmholtz free energy; difference between function and Euler 
# formulation
################################################################################
helmholtzFunction   = helmholtz(T,V_vap,n_vap)
helmholtzEuler      = pressure(T,V_vap,n_vap)*V_vap + dot(chemicalPotential(T,V_vap,n_vap),n_vap)

println("\nRelative error in Helmholtz free energy between function and Euler formulation:")
println((helmholtzFunction-helmholtzEuler)/helmholtzFunction)
println("\n")
