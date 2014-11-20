#############################################################
# Thermodynamic functions using ideal gas
#
# Functions for calculating thermodynamic potentials, 
# properties and useful derivatives using ideal gas EOS: 
#	- idealGasH: enthalpy for an ideal gas
#	- idealGasS: entropy for an ideal gas
#	- idealGasMu: chemical potential for an ideal gas
#	- idealGasA: Helmholtz free energy for an ideal gas
#
# Derivatives of Helmholtz free energy using ideal gas:
#	- Aig_T: 	$\pdc{A\ig}{T}{V,\vt{n}} = - S\ig$ 
#	- Aig_V: 	$\pdc{A\ig}{V}{T,\vt{n}} = - p\ig$
#	- Aig_n: 	$\pdc{A\ig}{\vt{n}}{T,V} = \mu\ig$
#	- Aig_TT:	$\pddc{A\ig}{T}{T}{V,\vt{n}}$ 
#	- Aig_TV:	$\pddc{A\ig}{T}{V}{\vt{n}}$ 
#	- Aig_Tn:	$\pddc{A\ig}{T}{\vt{n}}{V}$
#	- Aig_VV:	$\pddc{A\ig}{V}{V}{T,\vt{n}}$
#	- Aig_nV:	$\pddc{A\ig}{\vt{n}}{V}{T}$
#	- Aig_nn:	$\pddc{A\ig}{\vt{n}}{\vt{n}}{T,V}$
#
# Resulting Hessian matrix based on Helmholtz free energy
# for ideal gas with constant temperature $T$: 
# 	- idealGasHessian 
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-19 17:08:35
#############################################################

function idealGasH(T)
	# Calculate the enthalpy of an ideal gas using 
	# $\vt{h}\ig(T) = \Delta_{f} \vt{h}\ig(T_{\mathrm{ref}}) + \int_{T_{\mathrm{ref}}}^{T} \! \vt{c}_{p}(\tau) \, \mathrm{d}\tau$
	h = h_ref + c_p*(T-T_ref)
end

function idealGasS(T,V,n)
	# Calculate the entropy of an ideal gas using
	# $\vt{s}\ig(T,V,\vt{n}) = \vt{s}\ig(T_{\mathrm{ref}}, p_{\mathrm{ref}}) + \int_{T_{\mathrm{ref}}}^{T} \! \frac{\vt{c}_{p}(\tau)}{\tau} \, \mathrm{d}\tau - R \ln\left(\frac{\vt{n}RT}{Vp_{\mathrm{ref}}}\right)$
	s = s_ref + c_p*log(T/T_ref) - R*log((n*R*T)/(V*p_ref))
end

function idealGasMu(T,V,n)
	# Calculate the chemical potential of an ideal gas 
	# $\mathrm{\mu}\ig(T,V,\vt{n}) = \vt{h}\ig - T\vt{s}\ig$
	h 	= idealGasH(T)
	s 	= idealGasS(T,V,n)
	mu	= h - T*s 
end

function idealGasA(T,V,n)
	# Calculate the Helmholtz free energy of an ideal gas 
	# $A(T,V,\vt{n}) = -pV + \sum_{i=1}^n \! \mu_i N_i$
	mu	= idealGasMu(T,V,n)
	p	= idealGasEOS(T,V,n)
	A 	= -p*V + dot(mu,n)
end

#############################################################
# First derivatives of Helmholtz free energy
#############################################################

function Aig_T(T,V,n)
	# Calculate the derivative of Helmholtz free energy of 
	# an ideal gas with respect to temperature
	# $\pdc{A\ig}{T}{V,\vt{n}} = - S\ig$
	s 	= idealGasS(T,V,n) 
	A_T = -dot(n,s)
end

function Aig_V(T,V,n)
	# Calculate the derivative of Helmholtz free energy of 
	# an ideal gas with respect to volume
	# $\pdc{A\ig}{V}{T,\vt{n}} = - p\ig$
	A_T = -idealGasEOS(T,V,n)
end

function Aig_n(T,V,n)
	# Calculate the derivative of Helmholtz free energy of 
	# an ideal gas with respect to mole vector
	# $\pdc{A\ig}{\vt{n}}{T,V} = \vt{\mu}\ig$
	A_n = idealGasMu(T,V,n)
end

#############################################################
# Second derivatives of Helmholtz free energy
#############################################################

function Aig_TT(T,V,n)
	# Calculate the second derivative of Helmholtz free 
	# energy of an ideal gas with respect to temperature
	# $\pddc{A\ig}{T}{T}{V,\vt{n}} = -\vtt{n}\left(\frac{\vt{c}_p}{T} - \frac{R}{T}\vt{e}\right)$
	x = c_p/T
	y = (R/T)*ones(length(n))
	A_TT = -dot(n,(x-y))
end

function Aig_TV(T,V,n)
	# Calculate the second derivative of Helmholtz free energy 
	# of an ideal gas with respect to temperature and volume
	# $\pddc{A\ig}{T}{V}{\vt{n}} = -\frac{\vtt{n} R}{V}$
	A_TV = -sum(n)*(R/V)
end

function Aig_Tn(T,V,n)
	# Calculate the second derivative of Helmholtz free energy 
	# of an ideal gas with respect to temperature and mole vector
	# $\pddc{A\ig}{T}{\vt{n}}{V} = - \vt{s}\ig + R$
	A_TV = -idealGasS(T,V,n) + R*ones(length(n))
end

function Aig_VV(T,V,n)
	# Calculate the second derivative of Helmholtz free energy 
	# of an ideal gas with respect to volume
	# $\pddc{A\ig}{V}{V}{T,\vt{n}} = \frac{\vtt{n} RT}{V^2}$
	A_VV = sum(n)*(R*T/V^2)
end

function Aig_nV(T,V,n)
	# Calculate the second derivative of Helmholtz free energy 
	# of an ideal gas with respect to mole vector and volume
	# $\pddc{A\ig}{\vt{n}}{V}{T} = -\frac{\vt{e} RT}{V}$
	A_nV = -(ones(length(n)))*((R*T)/V)
end

function Aig_nn(T,V,n)
	# Calculate the second derivative of Helmholtz free energy 
	# of an ideal gas with respect to mole vector
	# $\pddc{A\ig}{\vt{n}}{\vtt{n}}{T,V} = -RT(\mathrm{diag}(\vt{n}))^{-1}$
	A_nn = -R*T*(inv(diagm(n)))
end

