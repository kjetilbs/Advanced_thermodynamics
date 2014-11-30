#############################################################
# Read experimental data
#
# Read experimental data from .csv-files, convert this into 
# SI-units and output data for plotting in (V,T)-diagram to 
# find suitable iteration area
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-19 14:05:42
#############################################################

dewPoints 		= readdlm("dewPoints.csv",',')
bubblePoints 	= readdlm("bubblePoints.csv",',')

# Pressure conversion: 6895 Pa = 1 psi
pressureConversion(p) = 6895*p

# Temperature conversion: Kelvin = (Fahrenheit - 32) * 5 / 9 + 273.15 
temperatureConversion(T) = ((T - 32*ones(length(T)))*(5/9))+273.15*ones(length(T))

# Density: given as g/cm^3 = 1000 kg/m^3
densityConversion(rho) = rho*1e3

# Feed composition of natural gas, ref. Gonzalez & Lee (1968)
n_feed 	= [0.016; 0.9450; 0.026; 0.0081; 0.0052];				# [-], mole fraction
Mm 		= [28.0134; 16.0425; 30.069; 44.0956; 58.1222]*1e-3; 	# [kg/mol], molar mass

# Convert the experimental data to SI units
dewPointT			= temperatureConversion(dewPoints[:,2]) 	# [K]
dewPointP			= pressureConversion(dewPoints[:,1])		# [Pa]
dewPointRho			= densityConversion(dewPoints[:,3])			# [kg/m3]

bubblePointT		= temperatureConversion(bubblePoints[:,2]) 	# [K]
bubblePointP		= pressureConversion(bubblePoints[:,1])		# [Pa]
bubblePointRho		= densityConversion(bubblePoints[:,3])		# [kg/m3]

# Calculating the weighted average molar mass
Mm_avg 				= dot(n_feed,Mm) 							# [kg/mol]

# Calculating the molar volume
dewPointV_m			= (1./dewPointRho)*Mm_avg 					# [m3/mol]
bubblePointV_m		= (1./bubblePointRho)*Mm_avg 				# [m3/mol]

# Output the data to .dat
writedlm("dewPoints_VT.dat",[dewPointV_m dewPointT], ' ')
writedlm("bubblePoints_VT.dat",[bubblePointV_m bubblePointT], ' ')

writedlm("dewPoints_VP.dat",[dewPointV_m dewPointP], ' ')
writedlm("bubblePoints_VP.dat",[bubblePointV_m bubblePointP], ' ')


