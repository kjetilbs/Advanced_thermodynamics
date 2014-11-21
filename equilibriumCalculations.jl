#############################################################
# Equilibrium calculation
#
# Equilibrium calculation for the phase equilibrium using the 
# Newton-Raphson method with thermodynamic data calculated 
# using the EOS module of choice.
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-21 15:43:32
#############################################################

module equilibriumCalculations

# Functions to be available in global scope once the module is
# imported/being used
export equilibriumCalculation

# Defining constants and reading component data. Needs to be 
# included within the module scope
include("defineConstants.jl")

# Using EOS modules
using idealGas
using redlichKwong

#############################################################
# Equilibrium calculation
#############################################################

function equilibriumCalculation(x, x_total, T)
	# Calculate the equilibrium state for a multiphase 
	# equilibrium problem with constant temperature $T$
	# given an initial guess and the constraints (in the
	# form of the total volume and total amount of substance)

	# Residual value; initialized to a large number
	residual 		= 1e10

	# Norm value; initialized to a large number
	norm 			= 1e10

	# Iteration counter; initial value
	iterationCount 	= 0

	# Maximum allowed iteration steps
	maxIterations 	= 1000

	# Convergence tolerance 
	convergenceTol 	= 1e-8

	# Convergence flag; 1 if converged, else 0
	hasConverged 	= false

	# Initializing iteration variables
	x_vap = x
	x_liq = x_total - x

	# Newton-Raphson iteration loop
	while !hasConverged && iterationCount < maxIterations
		# Printing iteration step
		println("Iteration nr.: "*string(iterationCount))

		# Unwrapping the state vector 
		V_vap 			= x_vap[1]
		n_vap 			= x_vap[2:end]

		V_liq 			= x_liq[1]
		n_liq 			= x_liq[2:end]

		# Calculating the current gradient vector for each phase 
		p_vap 			= pressure(T,V_vap,n_vap)
		mu_vap			= chemicalPotential(T,V_vap,n_vap)

		p_liq 			= pressure(T,V_liq,n_liq)
		mu_liq			= chemicalPotential(T,V_liq,n_liq)

		# Calculating the current Hessian matrix for each phase
		H_liq 			= hessian(T,V_vap,n_vap)
		H_vap 			= hessian(T,V_liq,n_liq)

		# Preparing the iteration
		currentGradient = [p_vap; mu_vap] - [p_liq; mu_liq]
		currentHessian 	= H_vap + H_liq

		# Applying the N-R iteration step
		deltaX 			= -currentHessian\currentGradient

		# Printing the current step
		println("Current step: "*string(deltaX))

		# Has converged? (i.e. sufficiently small change)
		if maximum(abs(deltaX)) < convergenceTol
			hasConverged = true
		end

		# Increment the iteration counter
		iterationCount += 1

		# Update the state according to the iteration step
		x_vap += deltaX
		x_liq -= deltaX

		# Printing the current state
		println("x_vap: "*string(x_vap))
		println("x_liq: "*string(x_liq))
		println("Deviation from total: "*string(x_total-(x_vap+x_liq)))
	end
end

#############################################################
# End module
#############################################################
end
