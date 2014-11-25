################################################################################
# Calculate numeric Hessian
#
# Julia function to calculate numeric Hessian to be used in the N-R-iteration in
# place of an analytic Hessian that (at the present time) does not work as it 
# should, presumably due to mistake(s).
#
# Author:     Kjetil Sonerud
# Updated:    2014-11-25 08:43:56
################################################################################

################################################################################
# Numeric Hessian calculation
################################################################################
function numericHessian(T,V_i,n_i)
    # Calculate a numeric Hessian from a central difference approximation, where
    # the volume and mole vector is perturbed in turn. Here, T is the current 
    # temperature, $V_i$ is the phase volume and $n_i$ is the phase mole vector. 

    # Calculating the current gradient vector
    p_i               = pressure(T,V_i,n_i)
    mu_i              = chemicalPotential(T,V_i,n_i) 
    gradient_i        = [p_i, mu_i]
        
    # Small volume perturbation
    delta_V_i         = 1e-11
    # Small mole number perturbation
    delta_n_i         = 1e-11
    
    # Initializing the numeric Hessian from perturbed gradient vectors
    numericH = zeros(1+length(n_i),1+length(n_i))

    # Preparing for central difference; volume
    forward_V_i       = V_i + 0.5*delta_V_i
    backward_V_i      = V_i - 0.5*delta_V_i
        
    forward_gradient_i    =     [
                                    -pressure(T,forward_V_i,n_i), 
                                    chemicalPotential(T,forward_V_i,n_i)
                                ]
    backward_gradient_i    =    [
                                    -pressure(T,backward_V_i,n_i), 
                                    chemicalPotential(T,backward_V_i,n_i)
                                ]
    
    # Result of volume perturbation
    numericH[1,:]       = centralDifference(forward_gradient_i, backward_gradient_i, delta_V_i)
    
    # Mole number perturbation
    for i in 1:length(n_i)
        # Preparing for central difference
        forward_n_i       = deepcopy(n_i)
        backward_n_i      = deepcopy(n_i)
    
        forward_n_i[i]    += 0.5*delta_n_i
        backward_n_i[i]   -= 0.5*delta_n_i
        
        forward_gradient_i    =     [
                                        -pressure(T,V_i,forward_n_i), 
                                        chemicalPotential(T,V_i,forward_n_i)
                                    ]
        backward_gradient_i    =    [
                                        -pressure(T,V_i,backward_n_i), 
                                        chemicalPotential(T,V_i,backward_n_i)
                                    ]
    
        # # Checking mole vectors
        # println("Before perturbation:   "*string(n_i))
        # println("Forward perturbation:  "*string(forward_n_i))
        # println("Backward perturbation: "*string(backward_n_i))
        # println("\n")

        # Adding to the numeric Hessian from perturbed gradient
        numericH[i+1,:] = centralDifference(forward_gradient_i, backward_gradient_i, delta_n_i)
    end
    # Output resulting numeric Hessian
    return numericH
end

################################################################################
# Central difference scheme
################################################################################
function centralDifference(f_forward, f_backward, h)
    # Using central difference scheme to estimate numeric
    # derivative using
    # $f'(x) \approx \frac{(f(x+0.5h) - f(x-0.5h))}{h}$
    dfun = (f_forward - f_backward)/h
end

