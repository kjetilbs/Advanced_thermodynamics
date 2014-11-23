#############################################################
# Legendre transform
#
# A short Julia snippet to illustrate the properties of the
# Legendre transform
#
# Author: 	Kjetil Sonerud
# Updated:	2014-11-18 10:30:52
#############################################################

#############################################################
# Example: f(x) = x^2 + 2
#############################################################
# Open file to write results to
filename = "plotting/legendre_original.dat"
f_handle = open(filename, "w")

# Define the original function
f(x) = x^2 + 2

# Define the domain
x_vec = linspace(-5,5,30)

# Calculate the function values of the original function
for i in x_vec 
	write(f_handle, string(i) * " " * string(f(i)) * "\n")
end

# Close file
close(f_handle)

#############################################################

# Define the derivative
df(x) = 2x

# Define the Legendre transform
legendre(theta) = -0.25theta^2 + 2

# Tangent line
y(theta, x) = theta*x + legendre(theta)

# Theta values
theta_vec = linspace(-10,10,30);

# Redefine the number of points in domain
x = linspace(-5,5,3)

# Open file to write results to
filename = "plotting/legendre_transform.dat"
f_handle = open(filename, "w")

# Tangent lines
for i in 1:length(theta_vec)
	# Calculating the tangent lines
	for j in x
		write(f_handle, string(j) * " " * string(y(theta_vec[i],j)) * "\n")
	end
	write(f_handle, "\n")
end

# Close file
close(f_handle)

#############################################################

# Open file to write results to
filename = "plotting/legendre_transform_fn.dat"
f_handle = open(filename, "w")

# Calculate the function values of the transformed function
for i in theta_vec
	write(f_handle, string(i) * " " * string(legendre(i)) * "\n")
end

# Close file
close(f_handle)

#############################################################

# Using a similar approach as above, we can see that using the
# original function to create tangent values will give the 
# Legendre-function in the theta-space

# Double Legendre 
legendre2(x) = x^2 + 2

# Tangent line
psi(theta, x) = x*theta + legendre2(x)

# Open file to write results to
filename = "plotting/legendre_original_tn.dat"
f_handle = open(filename, "w")

# Redefine the number of theta values
theta = linspace(-10,10,3);

# Tangent lines
for i in 1:length(x_vec)
	# Calculating the tangent lines
	for j in theta
		write(f_handle, string(j) * " " * string(psi(j,theta_vec[i])) * "\n")
	end
	write(f_handle, "\n")
end