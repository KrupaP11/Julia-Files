# SPS_Homework 5 Question 4
# Dr. Warren gave this code in python. Converting it to julia
#	in order to make it work with my code.

# implement this function below

function ridders(f, a, b, max_itrs, tol=1e-6)

	# Make sure that x0 < x1

	if a < b
		x₀ = a
		x₁ = b
		else
			x₀ = b
			x₁ = a
	end

	# Compute the function at the two endpoints and make sure 
	#	they bracket a root

	f₀ = f(x₀)
	f₁ = f(x₁)

	if f₀ == 0.0
		return x₀
	elseif f₁ == 0.0
		return x₁
	end

	# Throw an error if the root isn't bracketed
	if f₀*f₁ > 0.0
		error("Starting poitns for Ridder's rule must bracket the root")
	end

	# Main Ridders' method loop

	for i in 1:max_itrs
		# Find the midpoint of the interval

		x₂ = (x₀+x₁)/2.0

		f₂ = f(x₂)

		# Apply Ridders' formula to get x_3
		x₃ = x₂ + (x₁-x₂)*(f₂/f₀) / sqrt((f₂/f₀)^2 - (f₁/f₀))
		f₃ = f(x₃)
		#println(i, x₃)

		# Check to see if we can return: we took a very small step 
		#	or out lastest guess evaluates to 0
		if abs(x₃-x₂) < tol || f₃ == 0.0
			return x₃, i
		end

		# Bracket the root as tightly as possible
		if f₃*f₂ < 0.0
			# if this is the case then the root lies between x2 and x3, 
			#	so we set that as the interval for the next iteration

			x₀ = x₂; f₀ = f₂
			x₁ = x₃; f₁ = f₃
		else
			# f3 and f2 have the same sign so x3 adn x2 don't bracket the root.
			#	Check whether f0 or f1 is the other endpoint
			if f₃*f₁ < 0.0
				# x3 and x1 bracket the root, so replace x0 with x3
				x₀ = x₃; f₀ = f₃
			else
				# x3 and x0 bracket the root, so replace x1 with x3
				x₁ = x₃; f₁ = f₃
			end
		end
	end
	println("Exceeded iteration limit wihtout solution")
	return nothing

end