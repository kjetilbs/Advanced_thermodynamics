################################################################################
# Equilibrium calculation
#
# Equilibrium calculation for the phase equilibrium using the Newton-Raphson 
# method with thermodynamic data calculated using the EOS module of choice.
#
# Author:   Kjetil Sonerud
# Updated:  2014-11-30 18:24:49
################################################################################

module equilibriumCalculations

# Functions to be available in global scope 
export equilibriumCalculation, phaseEquilibrium

# Defining constants and reading component data
include("defineConstants.jl")

# Including numeric Hessian function
include("numericHessian.jl")

# Using EOS modules
# using idealGas
# using redlichKwong

################################################################################
# Inner loop: equilibrium calculation using the Newton-Raphson method
################################################################################

function equilibriumCalculation(x, x_total, T)
    # Calculate the equilibrium state for a multiphase 
    # equilibrium problem with constant temperature $T$
    # given an initial guess and the constraints (in the
    # form of the total volume and total amount of substance)

    # Iteration counter; initial value
    iterationCount  = 0

    # Maximum allowed iteration steps
    maxIterations   = 10000

    # Convergence tolerance 
    # convergenceTol[1]:    pressure tolerance
    # convergenceTol[2]:    chemical potential tolerance
    convergenceTol  = [5e-8, 1e-10]

    # Convergence flag; 1 if converged, else 0
    hasConverged    = false

    # Initializing iteration variables
    x_vap           = x
    x_liq           = x_total - x

    # Adjusting the volume guess
    x_liq[1] = 1.3*redlichKwong.redlichKwongB(x_liq[2:end])
    x_vap[1] = x_total[1] - x_liq[1]

    # Newton-Raphson iteration loop
    while !hasConverged && iterationCount < maxIterations
        # Printing iteration step
        if DEBUG
            println("\n")
            println("Iteration nr.: "*string(iterationCount))
        end

        # Unwrapping the state vector 
        V_vap           = x_vap[1]
        n_vap           = x_vap[2:end]

        V_liq           = x_liq[1]
        n_liq           = x_liq[2:end]

        # Calculating the current gradient vector for each phase 
        p_vap           = pressure(T,V_vap,n_vap)
        mu_vap          = chemicalPotential(T,V_vap,n_vap)

        p_liq           = pressure(T,V_liq,n_liq)
        mu_liq          = chemicalPotential(T,V_liq,n_liq)

        # Calculating the current Hessian matrix for each phase
        H_liq           = hessian(T,V_vap,n_vap)
        H_vap           = hessian(T,V_liq,n_liq)

        # Calculating the current Hessian matrix for each phase, using the 
        # numeric Hessian routine
        # H_liq           = numericHessian(T,V_vap,n_vap,"rk")
        # H_vap           = numericHessian(T,V_liq,n_liq,"rk")

        # Preparing the iteration
        currentGradient = [-p_vap, mu_vap] - [-p_liq, mu_liq]
        currentHessian  = H_vap + H_liq

        ########################################################################
        # Debugging gradient
        ########################################################################
        if DEBUG
            println("p_vap: "*string(p_vap))
            println("mu_vap: "*string(mu_vap))
            println("p_liq: "*string(p_liq))
            println("mu_liq: "*string(mu_liq))
            println("currentGradient: "*string(currentGradient))
        end
        ########################################################################

        # Applying the N-R iteration step
        deltaX          = -currentHessian\currentGradient

        # Copy current state to check the step length
        x_vap_try = copy(x_vap)
        x_liq_try = copy(x_liq)

        # Attempt to update the state according to the iteration step
        x_vap_try += deltaX
        x_liq_try -= deltaX

        # Counter for debugging
        stepReductionCounter = 0

        while ((checkSolution(x_vap_try, x_liq_try) == false) && (minimum(abs(deltaX)) > abs(1e-15)))
            # Take half the step length 
            deltaX = deltaX/2

            if DEBUG
                stepReductionCounter += 1
                println("Reducing step length... "*string(stepReductionCounter)". attempt")
            end

            # New copy of the current state
            x_vap_try = copy(x_vap)
            x_liq_try = copy(x_liq)

            # Update the copy of the current state according to time step
            x_vap_try += deltaX
            x_liq_try -= deltaX
        end

        # Updating the actual state after possible step size reduction 
        x_vap += deltaX
        x_liq -= deltaX

        # Has converged? (i.e. deviation from $\nabla A = 0$ sufficiently small)
        if norm(currentGradient[1]) < convergenceTol[1] && norm(currentGradient[2:end]) < convergenceTol[2]
            hasConverged = true
        end

        # Increment the iteration counter
        iterationCount += 1
    end

    # Check if the solution has converged
    if hasConverged == true
        # Equilibrium composition for the vapor phase and number of iterations
        return x_vap, iterationCount
    else
        # If no solution is found, display error to the user
        error("No solution is found for ("*string(x_total[1])*","*string(T)*")")
    end
end

################################################################################
# Equilibrium calculation (N-R-loop)
################################################################################

function phaseEquilibrium(x_guess, n_total, rangeT, rangeV)
    # Calculate the necessary phase equilibrium data for a multiphase equilibrium 
    # problem in the temperature range 
    # rangeT = [$T_{\mathrm{min}}$, $\ldots$ , $T_{\mathrm{max}}$] 
    # and volume range
    # rangeV = [$V_{\mathrm{min}}$, $\ldots$ , $V_{\mathrm{max}}$] 
    # given an initial guess and the constraints 
    
    # Dimensions of the iteration grid
    numVolumes          = length(rangeV)
    numTemperatures     = length(rangeT)

    # Initializing the solution vectors
    ansTemperature      = zeros(numVolumes, numTemperatures)
    ansPressure         = zeros(numVolumes, numTemperatures)
    ansVolumeTotal      = zeros(numVolumes, numTemperatures)
    ansVolumeVap        = zeros(numVolumes, numTemperatures)
    ansVolumeLiq        = zeros(numVolumes, numTemperatures)
    ansCompositionVap   = [zeros(length(n_total)) for x = 1:numVolumes, y = 1:numTemperatures]
    ansCompositionLiq   = [zeros(length(n_total)) for x = 1:numVolumes, y = 1:numTemperatures]
    ansEntropy          = zeros(numVolumes, numTemperatures)
    ansEnthalpy         = zeros(numVolumes, numTemperatures)
    ansIterations       = zeros(numVolumes, numTemperatures)

    # Debugging; storing the maximum error in $p$ and $\mu$
    maxNormP            = 0
    maxNormMu           = 0

    # Initializing temperature iteration array
    temperatureIterationRange   = [numTemperatures:-1:1]

    # Initializing volume iteration array
    volumeIterationRange        = [numVolumes:-1:1]

    # Initial guess for the iteration at ($T_{\mathrm{min}}$,$V_{\mathrm{min}}$), vapor phase
    #   - x[1]:         Vapor phase volume 
    #   - x[2:end]:     Vapor phase mole vector

    # Iterating on the volumes
    for volume in volumeIterationRange
        # To achieve the desired temperature range; flip the 
        # range at every iteration
        temperatureIterationRange = flipud(temperatureIterationRange)

        # Iterating on the temperature
        for temperature in temperatureIterationRange
            # Reset inner loop counter
            # innerLoopCounter    = 0

            # Preparing for the iteration
            V           = rangeV[volume]
            T           = rangeT[temperature]
            x_total     = [V, n_total]

            ####################################################################
            # Inner loop: 
            # Calculate the vapor-liquid equilibrium using Newton-Raphson method
            ####################################################################
            x_vap, innerLoopCounter = equilibriumCalculation(x_guess, x_total, T)
            x_liq                   = x_total - x_vap

            # Checking for trivial solution
            trivialCheck = x_vap[2:end]./x_total[2:end]

            # If trivial soultion found; warn the user and continue
            if round(maximum(trivialCheck),5) == round(minimum(trivialCheck),5)
                warn("Trivial solution found - will continue.")
                sleep(3)
            end

            println(
                        "\nIteration completed: ("*string(rangeT[temperature])*"K, "
                        *string(round(rangeV[volume],4))*" m3)"
                    )
            
            # Using this result as the initial guess for the next iteration
            x_guess     = copy(x_vap)

            # Unwrapping states
            V_vap       = x_vap[1]
            n_vap       = x_vap[2:end]
            V_liq       = x_liq[1]
            n_liq       = x_liq[2:end]

            # Checking the norm
            normP       =   norm(
                                    pressure(T,V_vap,n_vap) 
                                    - pressure(T,V_liq,n_liq)
                                )

            normMu      =   norm(
                                    chemicalPotential(T,V_vap,n_vap) 
                                    - chemicalPotential(T,V_liq,n_liq)
                                )

            # If norm is larger than current max; update
            maxNormP    = normP > maxNormP ? normP : maxNormP
            maxNormMu   = normMu > maxNormMu ? normMu : maxNormMu

            ####################################################################
            # Calculating the values of the remaining thermodynamic potentials 
            # at current $(T,V)$ and composition
            ####################################################################

            # Entropy
            S_vap       = entropy(T,V_vap,n_vap)
            S_liq       = entropy(T,V_liq,n_liq)
            ansEntropy[end-temperature+1, volume]               = S_vap + S_liq

            # Enthalpy
            H_vap       = enthalpy(T,V_vap,n_vap)
            H_liq       = enthalpy(T,V_liq,n_liq)
            ansEnthalpy[end-temperature+1, volume]              = H_vap + H_liq

            # Pressure
            ansPressure[end-temperature+1, volume]              = pressure(T,V_vap,n_vap)

            # Storing the variables at current $(T,V)$ and composition
            ansTemperature[end-temperature+1, volume]           = T
            ansVolumeTotal[end-temperature+1, volume]           = V
            ansVolumeVap[end-temperature+1, volume]             = V_vap
            ansVolumeLiq[end-temperature+1, volume]             = V_liq
            ansCompositionVap[end-temperature+1, volume][:]     = n_vap
            ansCompositionLiq[end-temperature+1, volume][:]     = n_liq
            ansIterations[end-temperature+1, volume]            = innerLoopCounter
        end
    end

    # Array of matrices containing the results
    ansArray = {ansTemperature, ansPressure, ansVolumeTotal, ansVolumeVap, ansVolumeLiq, 
                ansCompositionVap, ansCompositionLiq, ansEntropy, ansEnthalpy, ansIterations}

    ########################################
    # Display log data to the user
    ########################################
    println("\n\n")
    println("Max norm p: \n"*string(maxNormP))
    println("Max norm mu: \n"*string(maxNormMu))
    println("Max number of iterations: \n"*string(maximum(ansIterations)))
    println("Total number of iterations: \n"*string(sum(ansIterations)))
    println("Mean number of iterations: \n"*string(sum(ansIterations)/(numVolumes*numTemperatures)))
    println("\n\n")

    ########################################
    # Write log data to file
    ########################################
    fh_log   = open("results/log_"*string(strftime("%c", time()))*".txt", "w")

    # Writing data to log file
    write(fh_log, "Max norm p: " *string(maxNormP)*"\n")
    write(fh_log, "Max norm mu: " *string(maxNormMu)*"\n")
    write(fh_log, "Max number of iterations: " *string(maximum(ansIterations))*"\n")
    write(fh_log, "Total number of iterations: " *string(sum(ansIterations))*"\n")
    write(fh_log, "Mean number of iterations: " *string(sum(ansIterations)/(numVolumes*numTemperatures))*"\n")

    # Close files
    close(fh_log)
    ########################################

    return ansArray
end

################################################################################
# Auxiliary functions
################################################################################

function checkSolution(x_vap, x_liq)
    # Check that the solutions are physically meaningful:  
    #   - $V_{\mathrm{vap}} \:,\: V_{\mathrm{liq}} > B$
    #   - $\vt{n}_{\mathrm{vap}} \:,\: \vt{n}_{\mathrm{liq}} > 0$
    n_vap = x_vap[2:end]
    V_vap = x_vap[1]
    n_liq = x_liq[2:end]
    V_liq = x_liq[1]

    # Calculating $B$ for vapor and liquid phase
    B_vap = redlichKwong.redlichKwongB(n_vap)
    B_liq = redlichKwong.redlichKwongB(n_liq)

    # Assuring that the solutions are physically meaningful
    flag = (minimum(n_vap) > 0 && minimum(n_liq) > 0 && V_vap >= B_vap && V_liq >= B_liq)
end

################################################################################
# End module
################################################################################
end
