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
# Equilibrium calculation (N-R-loop)
#############################################################

function equilibriumCalculation(x, x_total, T)
	# Calculate the equilibrium state for a multiphase 
	# equilibrium problem with constant temperature $T$
	# given an initial guess and the constraints (in the
	# form of the total volume and total amount of substance)

	# Residual value; initialized to a large number
	residual 		= 1e10

	# Norm value; initialized to a large number
	# NB! Don't use norm, as it is a Julia function (!)
	myNorm 			= 1e10

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
		println("\n")

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

		# println(currentHessian)
		# hessianType = typeof(currentHessian)
		# println("Hessian type: "*string(hessianType))

		# println("Before...")
		# Printing the matrix norm
		matrixNorm = norm(currentHessian)
		# println("After!")

		println("Matrix norm: "*string(matrixNorm))
		println("\n")

		# Applying the N-R iteration step
		deltaX 			= -currentHessian\currentGradient

		# Printing the current step
		println("Current step: "*string(deltaX))
		println("\n")

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
		println("\n")
		println("\n")
		println("\n")
	end
end

#############################################################
# Equilibrium calculation (N-R-loop)
#############################################################

function phaseEquilibrium(x, x_total, rangeT, rangeV)
	# Calculate the necessary phase equilibrium data for a 
	# multiphase equilibrium problem in the temperature range 
	# rangeT = [$T_{\mathrm{min}}$; ... ; $T_{\mathrm{max}}$] 
	# and volume range
	# rangeV = [$V_{\mathrm{min}}$; ... ; $V_{\mathrm{max}}$] 
	# given an initial guess and the constraints 
	# (in the form of the total volume and total amount of substance)
	
	# Dimensions of the iteration grid
	numVolumes 		= length(rangeV)
	numTemperatures = length(rangeT)

	# Initializing the solution vectors
	ansPressure 	= zeros(numVolumes, numTemperatures)
	ansHelmholtz 	= zeros(numVolumes, numTemperatures)
	ansEntropy 		= zeros(numVolumes, numTemperatures)
	ansEnthalpy 	= zeros(numVolumes, numTemperatures)
	ansVolumeTotal 	= zeros(numVolumes, numTemperatures)
	ansVolumeVap 	= zeros(numVolumes, numTemperatures)
	ansVolumeLiq 	= zeros(numVolumes, numTemperatures)
	ansTemperature 	= zeros(numVolumes, numTemperatures)

	# push! vs preinitializing and counters?

	# Want to start from $T_{\mathrm{min}}$ at rangeV[1], then $T_{\mathrm{max}}$
	# at rangeV[2], the $T_{\mathrm{min}}$ at rangeV[3] and so 
	# on, in accordance with the proposed iteration scheme.
	# A suitable counter will take care of this.  
	volumeCounter = 0

	for volume in rangeV
		
		volumeCounter += 1
		if mod(volumeCounter,2) == 1
			rangeT = flipud(rangeT)
		end

		for temperature in rangeT
			# Do it!
		end
	end

end


#############################################################
# End module
#############################################################
end
