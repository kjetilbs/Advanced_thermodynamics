################################################################################
# Hessian testing
#
# Comparing the analytic Hessain and numeric Hessian (salculated using a simple 
# central difference scheme) for ideal gas, residual Redlich-Kwong and the actual 
# Redlich-Kwong equation of state. The goal is to check whether the analytic 
# expressian matches the numeric, with the latter taken to be correct.
#
# Author:     Kjetil Sonerud
# Updated:    2014-11-28 09:59:06
################################################################################

workspace()

########################################
# Importing EOS and numeric Hessian
########################################
# Using ideal gas module
include("idealGas.jl")
using idealGas

# Using Redlich-Kwong module
include("redlichKwong.jl")
using redlichKwong

# Numeric Hessian function
include("numericHessian.jl")

########################################
# Declaring variables
########################################
# Feed composition of natural gas, ref. Gonzalez & Lee (1968)
n_feed = [0.016 0.9450 0.026 0.0081 0.0052];

# Normalizing; converting to vector in Julia
n_feed = vec(n_feed/sum(n_feed))

# Distribution of components in the two phases
# n_vapor = [0.9; 0.8; 0.1; 0.01; 0.001].*n_feed
n_vap = [0.9; 0.8; 0.3; 0.05; 0.1].*n_feed

n_liq = n_feed - n_vap

# Calculating the hard sphere volume
V_liq = 1.3*redlichKwong.redlichKwongB(n_liq)
V_vap = 9e-4 - V_liq

# Temperature
T = 180

println("########################################")
println("# Ideal Hessian comparison")
println("########################################")

println("\nNumeric ideal Hessian:")
numericIdealHessian     = numericHessian(T,V_vap,n_vap,"ideal")
println(numericIdealHessian)
println("\nAnalytic ideal Hessian:")
analyticIdealHessian    = idealHessian(T,V_vap,n_vap)
println(analyticIdealHessian)
println("\nRelative error:")
println((numericIdealHessian-analyticIdealHessian)./numericIdealHessian)
println("\nNorm:")
println(norm((numericIdealHessian-analyticIdealHessian)./numericIdealHessian))

println("\n")

println("########################################")
println("# Residual Hessian comparison")
println("########################################")

println("\nNumeric residual Hessian:")
numericResidualHessian     = numericHessian(T,V_vap,n_vap,"res")
println(numericResidualHessian)
println("\nAnalytic residual Hessian:")
analyticResidualHessian    = redlichKwong.residualHessian(T,V_vap,n_vap)
println(analyticResidualHessian)
println("\nRelative error:")
relResError = (numericResidualHessian-analyticResidualHessian)./numericResidualHessian
println(round(relResError,6))
println("\nNorm:")
println(norm(relResError))
println("\nSymmetric relative error?")
println((relResError - relResError')./relResError)

println("\n")

println("########################################")
println("# Redlich-Kwong Hessian comparison")
println("########################################")

println("\nNumeric R-K Hessian:")
numericRKHessian     = numericHessian(T,V_vap,n_vap,"rk")
println(numericRKHessian)
println("\nAnalytic R-K Hessian:")
analyticRKHessian    = hessian(T,V_vap,n_vap)
println(analyticRKHessian)
println("\nRelative error:")
println(round((numericRKHessian-analyticRKHessian)./numericRKHessian,6))
println("\nNorm:")
println(norm((numericRKHessian-analyticRKHessian)./numericRKHessian))



