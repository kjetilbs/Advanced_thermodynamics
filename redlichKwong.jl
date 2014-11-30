################################################################################
# Redlich-Kwong module
#
# Redlich-Kwong equation of state:
#   - redlichKwongA:    Parameter A used in the R-K EOS
#   - redlichKwongB:    Parameter B used in the R-K EOS
#   - pressure:         Pressure from the R-K EOS
#
# Functions for calculating thermodynamic potentials, 
# properties and useful derivatives using Redlich-Kwong EOS: 
#
# Residual functions using R-K EOS: 
#   - residualHelmholtz:            Residual Helmholtz free energy
#   - residualEntropy:              Residual entropy
#   - residualChemicalPotential:    Residual chemical potential
#   - residualPressure:             Residual pressure, $p\rrk \eqdef p\rk - p\ig$
# 
# Thermodynamic potentials using Redlich-Kwong EOS, found by combining the 
# residual potential and the ideal gas potential:
#   - helmholtz:            $A\rk(T,V,\vt{n}) = A\rrk(T,V,\vt{n}) + A\ig(T,V,\vt{n})$
#   - entropy:              $S\rk(T,V,\vt{n}) = S\rrk(T,V,\vt{n}) + S\ig(T,V,\vt{n})$
#   - chemicalPotential:    $\mu\rk(T,V,\vt{n}) = \mu\rrk(T,V,\vt{n}) + \mu\ig(T,V,\vt{n})$
#   - enthalpy:             $H\rk(T,V,\vt{n}) = A\rk(T,V,\vt{n}) + TS\rk(T,V,\vt{n}) + p\rk(T,V,\vt{n})V$
#
# First derivatives of residual Helmholtz free energy:
#   - residualA_V:      $\pdc{A\rrk}{V}{T,\vt{n}} \eqdef - (p\rk - p\ig)$ 
#
# Second derivatives of residual Helmholtz free energy:
#   - residualA_VV:     $\pddc{A\rrk}{V}{V}{T,\vt{n}}$
#   - residualA_nV:     $\pddc{A\rrk}{\vt{n}}{V}{T}$
#   - residualA_nn:     $\pddc{A\rrk}{\vt{n}}{\vtt{n}}{T,V}$ 
# 
# Resulting Hessian matrix based on residual Helmholtz free energy using R-K, 
# and the R-K Hessian as the sum of the ideal gas Hessian and the residual
# Hessian; all with constant temperature $T$: 
#   - residualHessian:  $\mt{H}\rrk$
#   - hessian:          $\mt{H}\rk$
#
# Author:   Kjetil Sonerud
# Updated:  2014-11-30 17:43:44
################################################################################

module redlichKwong

# Functions to be available in global scope 
export pressure, entropy, enthalpy, chemicalPotential, helmholtz, hessian

# Defining constants and reading component data
include("defineConstants.jl")

# Using the ideal gas module
# include("idealGas.jl")
using idealGas

################################################################################
# Redlich-Kwong EOS
################################################################################
function redlichKwongA(T,n)
    # Calculate the parameter A for the Redlich-Kwong EOS
    # $A(T,\vt{n}) = \frac{1}{T^{1/2}}\vtt{n}(\vt{a}\:\vtt{a})\vt{n}$
    # where $a_i$ is defined in defineConstants.jl as 
    # $a_i = \frac{1}{9(2^{1/3}-1)}\frac{R^2T_{c,i}^{5/2}}{p_{c,i}} = \Omega_a \frac{R^2T_{c,i}^{5/2}}{p_{c,i}}$
    A = dot(n,((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T))
end

function redlichKwongB(n)
    # Calculate the parameter B for the Redlich-Kwong EOS
    # $B(\vt{n}) = \vtt{b}\:\vt{n}$
    # where $b_i$ is defined in defineConstants.jl as
    # $b_i = \frac{2^{1/3}-1}{3}\frac{RT_{c,i}}{p_{c,i}} = \Omega_b \frac{RT_{c,i}}{p_{c,i}}$ 
    B = dot(b_RK,n)
end

function pressure(T,V,n)
    # Calculate the pressure from the Redlich-Kwong EOS:
    # $p\rk = (T,V,\vt{n}) = \frac{NRT}{V-B} - \frac{A}{V(V+B)}$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    p = (sum(n)*R*T)/(V-B) - (A/(V*(V+B)))
end

################################################################################
# Residual thermodynamic potentials using Redlich-Kwong EOS
################################################################################
function residualHelmholtz(T,V,n)
    # Calculate the residual Helmholtz free energy using the R-K EOS
    # $\text{See \cref{eq:Helmholtz_residual_derivation_5}}$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    resA = sum(n)*R*T*log(V/(V-B)) + (A/B)*log(V/(V+B))
end

function residualEntropy(T,V,n)
    # Calculate the residual entropy using the R-K EOS
    # $\text{See \cref{eq:S_rrk}}$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    resS = -sum(n)*R*log(V/(V-B)) + (A/(2*B*T))*log(V/(V+B))
end

function residualChemicalPotential(T,V,n)
    # Calculate the residual chemical potential using the R-K EOS
    # $\text{See \cref{eq:mu_rrk}}$

    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    # Calculating $\Adn$ and $\Bdn$ parameters
    A_n = 2*((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)
    B_n = b_RK

    resMu =     (
                    ones(n)*R*T*log(V/(V-B)) + sum(n)*R*T*(B_n/(V-B)) 
                    + ((A_n*B - B_n*A)/(B^2))*log(V/(V+B)) - (A/B)*(B_n/(V+B))
                )
end

function residualPressure(T,V,n)
    # Calculate the residual pressure from the Redlich-Kwong EOS and ideal gas:
    # $p\rrk \eqdef p\rk - p\ig$
    resP = pressure(T,V,n) - idealPressure(T,V,n)
end

################################################################################
# Thermodynamic potentials using Redlich-Kwong EOS:
# combining the residual potential and the ideal gas potential
################################################################################
function helmholtz(T,V,n)
    # Calculate the Helmholtz free energy using the R-K EOS
    # $A\rk(T,V,\vt{n}) = A\rrk(T,V,\vt{n}) + A\ig(T,V,\vt{n})$
    rkA = residualHelmholtz(T,V,n) + idealHelmholtz(T,V,n)
end

function entropy(T,V,n)
    # Calculate the entropy using the R-K EOS
    # $S\rk(T,V,\vt{n}) = S\rrk(T,V,\vt{n}) + S\ig(T,V,\vt{n})$
    rkS = residualEntropy(T,V,n) + dot(n,idealEntropy(T,V,n))
end

function chemicalPotential(T,V,n)
    # Calculate the entropy using the R-K EOS
    # $\mu\rk(T,V,\vt{n}) = \mu\rrk(T,V,\vt{n}) + \mu\ig(T,V,\vt{n})$
    rkMu = residualChemicalPotential(T,V,n) + idealChemicalPotential(T,V,n)
end

function enthalpy(T,V,n)
    # Calculate the enthalpy using the R-K EOS from
    # $H\rk(T,V,\vt{n}) = A\rk(T,V,\vt{n}) + TS\rk(T,V,\vt{n}) + p\rk(T,V,\vt{n})V$
    rkH = helmholtz(T,V,n) + T*entropy(T,V,n) + pressure(T,V,n)*V
end

################################################################################
# First derivatives of residual Helmholtz free energy using 
# Redlich-Kwong EOS
################################################################################
function residualA_V(T,V,n)
    # Calculate the first derivative of $A\rrk$ with respect to $V$
    # $\pdc{A\rrk}{V}{T,\vt{n}} \eqdef - (p\rk - p\ig)$
    resA_V = -(pressure(T,V,b) - idealP(T,V,n))
end

################################################################################
# Second derivatives of residual Helmholtz free energy using 
# Redlich-Kwong EOS
################################################################################
function residualA_VV(T,V,n)
    # Calculate the second derivative of $A\rrk$ with respect to $V$
    # $\text{See \cref{eq:dA_rrk_dVdV}}$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    resA_VV = sum(n)*R*T*((2*V*B - B^2)/(V^2*(V-B)^2)) - A*((2*V + B)/(V^2*(V+B)^2))
end

function residualA_nV(T,V,n)
    # Calculate the second derivative of $A\rrk$ with respect to $\vt{n}$ and $V$
    # $\text{See \cref{eq:dA_rrk_dndV}}$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    # Calculating $\Adn$ and $\Bdn$ parameters
    A_n = 2*((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)
    B_n = b_RK

    resA_nV =   (
                    ones(length(n))*((R*T)/V) - ((ones(length(n))*R*T*(V-B) + B_n*sum(n)*R*T)/(V-B)^2)
                    + (A_n*B*V + B^2*A_n - A*B*B_n)/(B*V*(V+B)^2)
                )
end

function residualA_nn(T,V,n)
    # Calculate the second derivative of $A\rrk$ with respect $\vt{n}$ 
    # $\text{See \cref{eq:dA_rrk_dnn}}$

    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    # Calculating $\Adn$ and $\Bdn$ parameters
    A_n = 2*((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)
    B_n = b_RK

    # Calculating $\Adnn$ parameter
    A_nn = 2*((sqrt(a_RK)*sqrt(a_RK)'))/sqrt(T)

    # Temporary variables for the terms
    term1   = ((R*T)/(V-B))*((ones(n)*B_n') + (B_n*(ones(n)')))

    term2   = sum(n)*R*T*((B_n*B_n')/((V-B)^2))

    term3   = ((B^3*A_nn-B^2*((A_n*B_n')+(B_n*A_n'))+2*A*B*(B_n*B_n'))/(B^4))*log(V/(V+B))

    term4_1 = 2*A*V*(B_n*B_n')-V*B*((A_n*B_n')+(B_n*A_n'))+3*A*B*(B_n*B_n')
    
    term4_2 = - B^2*((A_n*B_n')+(B_n*A_n'))
    
    term4   = (term4_1 + term4_2)/(B^2*(V+B)^2)

    # Result; summation of temporary terms
    resA_nn = term1 + term2 + term3 + term4
end

################################################################################
# Hessian matrices
################################################################################

function residualHessian(T,V,n)
    # Calculate the Hessian matrix based on residual Helmholtz 
    # free energy using R-K EOS with constant temperature $T$:
    A_VV    = residualA_VV(T,V,n)
    A_nV    = residualA_nV(T,V,n) 
    A_nn    = residualA_nn(T,V,n)

    H = [A_VV transpose(A_nV); A_nV A_nn]
end

function hessian(T,V,n)
    # Calculate the Hessian matrix based on Helmholtz free energy 
    # using R-K EOS with constant temperature $T$:
    H = residualHessian(T,V,n) + idealHessian(T,V,n)
end

################################################################################
# End module
################################################################################
end