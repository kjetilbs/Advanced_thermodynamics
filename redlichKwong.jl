################################################################################
# Thermodynamic functions using Redlich-Kwong
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
#   - residualA_TT:     $\pddc{A\rrk}{T}{T}{V,\vt{n}}$
#   - residualA_TV:     $\pddc{A\rrk}{T}{V}{V,\vt{n}}$
#   - residualA_Tn:     $\pddc{A\rrk}{T}{\vt{n}}{V}$
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
# Updated:  2014-11-20 11:28:33
################################################################################

module redlichKwong

# Functions to be available in global scope once the module is
# imported/being used
export pressure, entropy, enthalpy, chemicalPotential, helmholtz, hessian

# Defining constants and reading component data. Needs to be 
# included within the module scope
include("defineConstants.jl")

# Using the ideal gas module
# include("idealGas.jl")
using idealGas
# println("Contents of ideal gas module:")
# whos(idealGas)
# println("\n")

################################################################################
# Redlich-Kwong EOS
################################################################################
function redlichKwongA(T,n)
    # Calculate the parameter A for the Redlich-Kwong EOS
    # $A(T,\vt{n}) = \frac{1}{T^{1/2}}\vtt{n}(\vt{a}\:\vtt{a})\vt{n}$
    # where a is defined in defineConstants.jl as 
    # $a_i = \frac{1}{9(2^{1/3}-1)}\frac{R^2T_{c,i}^{5/2}}{p_{c,i}} = \Omega_a \frac{R^2T_{c,i}^{5/2}}{p_{c,i}}$
    A = dot(n,((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T))
end

function redlichKwongB(n)
    # Calculate the parameter B for the Redlich-Kwong EOS
    # $B(\vt{n}) = \vtt{b}\:\vt{n}$
    # where b is defined in defineConstants.jl as
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
    # Calculate the residual Helmholtz free energy using the
    # R-K EOS
    # $A\rrk(T,V,\vt{n}) = NRT\ln \left(\frac{V}{V-B}\right)+\frac{A}{B}\ln \left(\frac{V}{V+B}\right)$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    resA = sum(n)*R*T*log(V/(V-B)) + (A/B)*log(V/(V+B))
end

function residualEntropy(T,V,n)
    # Calculate the residual entropy using the R-K EOS
    # $S\rrk(T,V,\vt{n}) = -\pdc{A\rrk}{T}{V,\vt{n}} = -NR\ln \left(\frac{V}{V-B}\right)+\frac{A}{2BT}\ln \left(\frac{V}{V+B}\right)$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    resS = -sum(n)*R*log(V/(V-B)) - (A/(2*B*T))*log(V/(V+B))
end

function residualChemicalPotential(T,V,n)
    # Calculate the residual chemical potential using the 
    # R-K EOS
    # $\mu\rrk\tvn = \pdc{A\rrk}{\vt{n}}{\tvni} = \vt{e}RT\ln \left(\frac{V}{V-B}\right) + NRT \left(\frac{B_{\vt{n}}}{V-B}\right) + \left(\frac{\Adn B - \Bdn A}{B^2}\right)\ln \left(\frac{V}{V+B}\right) - \frac{A}{B}\frac{\Bdn}{V+B}$
    # where
    # $\Adn &\meqdef \pdc{A}{\vt{n}}{\tvni} = \frac{2}{T^{1/2}}(\vt{a}\:\vtt{a})\vt{n}$
    # $\Bdn &\meqdef \pdc{B}{\vt{n}}{\tvni} = \vt{b}$

    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    # Calculating $\Adn$ and $\Bdn$ parameters
    A_n = 2*((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)
    B_n = b_RK

    resMu = ones(length(n))*R*T*log(V/(V-B)) + sum(n)*R*T*(B_n/(V-B)) + ((A_n*B - B_n*A)/(B^2))*log(V/(V+B)) - (A/B)*(B_n/(V+B))
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
    # Calculate the first derivative of residual Helmholtz 
    # free energy with respect to volume
    # $\pdc{A\rrk}{V}{T,\vt{n}} \eqdef - (p\rk - p\ig)$
    resA_V = -(pressure(T,V,b) - idealP(T,V,n))
end

################################################################################
# Second derivatives of residual Helmholtz free energy using 
# Redlich-Kwong EOS
################################################################################
function residualA_TT(T,V,n)
    # Calculate the second derivative of residual Helmholtz 
    # free energy with respect to temperature
    # $\pddc{A\rrk}{T}{T}{V,\vt{n}} = \frac{3A}{4BT^2}\ln \left(\frac{V}{V+B}\right)$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    resA_TT = ((3*A)/(4*B*T^2))*log(V/(V+B))
end

function residualA_TV(T,V,n)
    # Calculate the second derivative of residual Helmholtz 
    # free energy with respect to temperature and volume
    # $\pddc{A\rrk}{T}{V}{V,\vt{n}} = \frac{A_T}{V(V+B)} -\frac{NRB}{V(V-B)}$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    # Calculating $A_T$ parameter
    A_T = -A/2*T

    resA_TV = (A_T)/(V*(V+B)) - (sum(n)*R*B)/(V*(V-B))
end

function residualA_Tn(T,V,n)
    # Calculate the second derivative of residual Helmholtz 
    # free energy with respect to temperature and mole vector
    # $\pddc{A\rrk}{T}{\vt{n}}{V} = \vt{e}R\ln \left(\frac{V}{V-B}\right) + NR \left(\frac{B_{\vt{n}}}{V-B}\right) + \left(\frac{-\frac{\Adn B}{2T} - \Bdn A_T}{B^2}\right)\ln \left(\frac{V}{V+B}\right) - \frac{A_T}{B}\frac{\Bdn}{V+B}$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    # Calculating $\Adn$ and $\Bdn$ parameters
    A_n = 2*((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)
    B_n = b_RK

    # Calculating $A_T$ parameter
    A_T = -A/2*T

    resA_Tn =   ones(length(n))*R*log(V/(V-B)) + sum(n)*R*(B_n/(V-B))+((-((A_n*B)/2*T) 
                - B_n*A_T)/(B^2))*log(V/(V+B))-((A_T*B_n)/(B*(V+B)))
end

function residualA_VV(T,V,n)
    # Calculate the second derivative of residual Helmholtz 
    # free energy with respect to volume
    # $\pddc{A\rrk}{V}{V}{T,\vt{n}} = NRT \left[-\frac{1}{V^2} + \frac{1}{(V-B)^2}\right]+\frac{A}{B}\left[-\frac{1}{V^2} + \frac{1}{(V+B)^2}\right]$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    resA_VV = sum(n)*R*T*((2*V*B - B^2)/(V^2*(V-B)^2)) - A*((2*V + B)/(V^2*(V+B)^2))
end

function residualA_nV(T,V,n)
    # Calculate the second derivative of residual Helmholtz 
    # free energy with respect to mole vector and volume
    # $\pddc{A\rrk}{\vt{n}}{V}{T} = \frac{\vt{e}RT}{V} - \frac{\vt{e}RT(V-B) + \Bdn NRT}{(V-B)^2} + \frac{(\Adn BV + B^2\Adn - AB\Bdn)}{BV(V+B)^2}$
    
    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    # Calculating $\Adn$ and $\Bdn$ parameters
    A_n = 2*((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)
    B_n = b_RK

    # Calculating $A_T$ parameter
    A_T = -A/2*T

    resA_nV =   (
                    ones(length(n))*((R*T)/V) - ((ones(length(n))*R*T*(V-B) + B_n*sum(n)*R*T)/(V-B)^2)
                    + (A_n*B*V + B^2*A_n - A*B*B_n)/(B*V*(V+B)^2)
                )
    # + ((A_n*B*V - A*B_n*V)/(B^2*V^2)) - ((A_n*B*(V+B) - A*(V*B_n+2*B*B_n))/(B^2*(V+B)^2))
end

function residualA_nn(T,V,n)
    # Calculate the second derivative of residual Helmholtz 
    # free energy with respect to mole vector 
    # $\pddc{A\rrk}{\vt{n}}{\vtt{n}}{T,V} = \frac{RT}{(V-B)}(\vt{e}\Bdn\trans + \Bdn\vtt{n}) + NRT\left(\frac{\Bdn\Bdn\trans}{(V-B)^2}\right) \nonumber \left(\frac{\:\Adnn B^3 - (\Adn\Bdn\trans + \Bdn\Adn\trans)B^2 + 2B\Bdn\Bdn\trans A}{B^4}\right)\ln \left(\frac{V}{V+B}\right)$
    #       $\left(\frac{(\Adn B - \Bdn A)}{B^2}\right)\left(\frac{\Bdn\trans}{(V+B)}\right) - \left(\frac{\Bdn}{(V+B)}\right)\left(\frac{(\Adn\trans B - \Bdn\trans A)}{B^2}\right) \left(\frac{A}{B}\right)\left(\frac{\Bdn\Bdn\trans}{(V+B)^2}\right)$
    # $\text{See \cref{eq:dA_rrk_dnn}}$

    # Calculating $A$ and $B$ parameters
    A = redlichKwongA(T,n)
    B = redlichKwongB(n)

    # Calculating $\Adn$ and $\Bdn$ parameters
    A_n = 2*((sqrt(a_RK)*sqrt(a_RK)')*n)/sqrt(T)
    B_n = b_RK

    # Calculating $\Adnn$ parameter
    A_nn = 2*((sqrt(a_RK)*sqrt(a_RK)'))/sqrt(T)

    # Initializing 
    # resA_nn = zeros(n*n')

    # for i = 1:length(n)
    #     for j = 1:length(n)
    #         resA_nn[i,j] =   (
    #                             (R*T*(b_RK[i]+b_RK[j]))/(V-B) + (sum(n)*R*T*b_RK[i]*b_RK[j])/(V-B)^2
    #                             + (2*sqrt(a_RK[i]/sqrt(T)*a_RK[j]/sqrt(T))/B
    #                             - (b_RK[j]*2*sqrt(a_RK[i]*A/sqrt(T))
    #                             + b_RK[i]*2*sqrt(a_RK[j]*A/sqrt(T)))/B^2 
    #                             + 2*A*b_RK[i]*b_RK[j]/B^3)*log(V/(V-B))
    #                             - (2*sqrt(a_RK[i]*A/sqrt(T))*b_RK[j]
    #                             + 2*sqrt(a_RK[j]*A/sqrt(T))*b_RK[i])/(B*(V+B))
    #                             + (2*A*b_RK[i]*b_RK[j]*V)/(B^2*(V+B)^2)
    #                             + (3*A*b_RK[i]*b_RK[j])/(B*(V+B)^2)
    #                         )
    #     end
    # end

    resA_nn =   (
                    ((R*T)/(V-B))*(ones(length(n))*B_n') + (B_n*ones(length(n))') 
                    + sum(n)*R*T*((B_n*B_n')/((V-B)^2))
                    + ((A_nn*B^3-((A_n*B_n')+(B_n*A_n'))*B^2+2*B*(B_n*B_n')*A)/(B^4))*log(V/(V+B))
                    - ((A_n*B - B_n*A)/(B^2))*((B_n')/(V+B)) 
                    - ((B_n)/(V+B))*(((A_n'*B) - (B_n'*A))/(B^2))
                    + (A/B)*((B_n*B_n')/(V+B)^2)
                )
    
    # # Defining terms
    # term1   = ((R*T)/(V-B))*(ones(length(n))*B_n') + (B_n*ones(length(n))')
    # term2   = sum(n)*R*T*((B_n*B_n')/((V-B)^2))
    # term3   = ((B^2*A_nn-B*((A_n*B_n')+(B_n*A_n'))+2*A*(B_n*B_n'))/(B^3))*log(V/(V+B))
    # term4_1 = 2*A*V*(B_n*B_n')-V*B*((A_n*B_n')+(B_n*A_n'))+A*B*(B_n*B_n')
    # term4_2 = - B^2*((A_n*B_n')+(B_n*A_n'))
    # term4   = (term4_1 + term4_2)/(B^2*(V+B)^2)

    # resA_nn = term1 + term2 + term3 + term4 

    # resA_nn =   (
    #                 ((R*T)/(V-B))*(ones(length(n))*B_n') + (B_n*ones(length(n))')
    #               + sum(n)*R*T*((B_n*B_n')/((V-B)^2))
    #               + ((B^2*A_nn - B*((A_n*B_n')+(B_n*A_n')) + 2*A*B*(B_n*B_n'))/(B^3))
    #               * log(V/(V+B))
    #               + (2*A*V*(B_n*B_n')-V*B*((A_n*B_n')+(B_n*A_n'))+A*B*(B_n*B_n')
    #               - B^2*((A_n*B_n')+(B_n*A_n')))/(B^2*(V+B)^2)
    #             )
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