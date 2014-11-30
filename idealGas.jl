################################################################################
# Ideal gas module
#
# Ideal gas equation of state:
#   - idealPressure:            pressure from the ideal gas law
#
# Functions for calculating thermodynamic potentials, 
# properties and useful derivatives using ideal gas EOS: 
#   - idealH:                   enthalpy for an ideal gas
#   - idealEntropy:             entropy for an ideal gas
#   - idealChemicalPotential:   chemical potential for an ideal gas
#   - idealHelmholtz:           Helmholtz free energy for an ideal gas
#
# Derivatives of Helmholtz free energy using ideal gas:
#   - Aig_T:    $\pdc{A\ig}{T}{V,\vt{n}} = - S\ig$ 
#   - Aig_V:    $\pdc{A\ig}{V}{T,\vt{n}} = - p\ig$
#   - Aig_n:    $\pdc{A\ig}{\vt{n}}{T,V} = \mu\ig$
#   - Aig_VV:   $\pddc{A\ig}{V}{V}{T,\vt{n}}$
#   - Aig_nV:   $\pddc{A\ig}{\vt{n}}{V}{T}$
#   - Aig_nn:   $\pddc{A\ig}{\vt{n}}{\vt{n}}{T,V}$
#
# Resulting Hessian matrix based on Helmholtz free energy
# for ideal gas with constant temperature $T$: 
#   - idealHessian 
#
# Author:   Kjetil Sonerud
# Updated:  2014-11-30 17:28:43
################################################################################

module idealGas

# Functions to be available in global scope
export  idealPressure, idealEntropy, idealChemicalPotential, idealHelmholtz, 
        idealHessian

# Defining constants and reading component data
include("defineConstants.jl")

################################################################################
# Ideal gas EOS
################################################################################
function idealPressure(T,V,n)
    # Calculate the pressure from the ideal gas EOS:
    # $p\ig(T,V,\vt{n}) = \frac{NRT}{V}$
    p = (sum(n)*R*T)/V
end

################################################################################
# Thermodynamic potentials using ideal gas
################################################################################
function idealH(T)
    # Calculate the enthalpy of an ideal gas using 
    # $\vt{h}\ig(T) = \Delta_{f} \vt{h}\ig(T_{\mathrm{ref}}) + \int_{T_{\mathrm{ref}}}^{T} \! \vt{c}_{p}(\tau) \, \mathrm{d}\tau$
    h = h_ref + c_p*(T-T_ref)
end

function idealEntropy(T,V,n)
    # Calculate the entropy of an ideal gas using
    # $\vt{s}\ig(T,V,\vt{n}) = \vt{s}\ig(T_{\mathrm{ref}}, p_{\mathrm{ref}}) + \int_{T_{\mathrm{ref}}}^{T} \! \frac{\vt{c}_{p}(\tau)}{\tau} \, \mathrm{d}\tau - R \ln\left(\frac{\vt{n}RT}{Vp_{\mathrm{ref}}}\right)$
    s = s_ref + c_p*log(T/T_ref) - R*log((n*R*T)/(V*p_ref))
end

function idealChemicalPotential(T,V,n)
    # Calculate the chemical potential of an ideal gas 
    # $\mathrm{\mu}\ig(T,V,\vt{n}) = \vt{h}\ig - T\vt{s}\ig$
    h   = idealH(T)
    s   = idealEntropy(T,V,n)
    mu  = h - T*s 
end

function idealHelmholtz(T,V,n)
    # Calculate the Helmholtz free energy of an ideal gas 
    # $A(T,V,\vt{n}) = -p\ig V  + \vt{n}\trans \vt{\mu}\ig$
    mu  = idealChemicalPotential(T,V,n)
    p   = idealPressure(T,V,n)
    A   = -p*V + dot(mu,n)
end

################################################################################
# First derivatives of Helmholtz free energy
################################################################################

function Aig_T(T,V,n)
    # Calculate the derivative of Helmholtz free energy of 
    # an ideal gas with respect to temperature
    # $\pdc{A\ig}{T}{V,\vt{n}} = - S\ig$
    s   = idealEntropy(T,V,n) 
    A_T = -dot(n,s)
end

function Aig_V(T,V,n)
    # Calculate the derivative of Helmholtz free energy of 
    # an ideal gas with respect to volume
    # $\pdc{A\ig}{V}{T,\vt{n}} = - p\ig$
    A_T = -idealPressure(T,V,n)
end

function Aig_n(T,V,n)
    # Calculate the derivative of Helmholtz free energy of 
    # an ideal gas with respect to mole vector
    # $\pdc{A\ig}{\vt{n}}{T,V} = \vt{\mu}\ig$
    A_n = idealChemicalPotential(T,V,n)
end

################################################################################
# Second derivatives of Helmholtz free energy
################################################################################

function Aig_VV(T,V,n)
    # Calculate the second derivative of Helmholtz free energy 
    # of an ideal gas with respect to volume
    # $\text{See \cref{eq:dA_ig_dVdV}}$
    A_VV = sum(n)*(R*T/V^2)
end

function Aig_nV(T,V,n)
    # Calculate the second derivative of Helmholtz free energy 
    # of an ideal gas with respect to mole vector and volume
    # $\text{See \cref{eq:dA_ig_dndV}}$
    A_nV = -(ones(n))*((R*T)/V)
end

function Aig_nn(T,V,n)
    # Calculate the second derivative of Helmholtz free energy 
    # of an ideal gas with respect to mole vector
    # $\text{See \cref{eq:dA_ig_dndn}}$
    A_nn = R*T*(diagm(1./n))
end

function idealHessian(T,V,n)
    # Calculate the Hessian matrix based on Helmholtz free energy
    # for ideal gas with constant temperature $T$:
    A_VV    = Aig_VV(T,V,n)
    A_nV    = Aig_nV(T,V,n) 
    A_nn    = Aig_nn(T,V,n)

    H = [A_VV transpose(A_nV); A_nV A_nn]
end

################################################################################
# End module
################################################################################
end