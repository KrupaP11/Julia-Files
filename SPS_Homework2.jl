# SPS 5100 Research Techniques - Homework 2
using Plots

function main()
	# Creating output folder
	# specify the directory name
	dir = "./output"

	# check if it exists, if not, create it
	if !isdir(dir)
	    mkpath(dir)
	    println("Directory created: $dir")
	end

	println("\nSolving question 1\n")
	derivative()

	println("\nSolving question 2\n")
	sick()

	println("\nSolving question 3\n")
	different_derivatives()

	println("\nQuestion 4 is defined on paper")

	println("\nSolving Question 5\n")
	partial_derivatives()

end

# Derivative function
function derivative()
	
	print("Solving part (a)\n")

	function f(x)
		return sin(π*x)
	end
	
	x = 0:2
	h = 0.01

	error = Float64[]
	for i in x
		cdf = ( f(i+h) - f(i-h) )/ (2*h)
		
		abs_error = cdf - (π*cos(π*i))
		push!(error, abs_error)
	end

	p = plot(x, error, label="Absolute error", lc=:black, lw=2)
	
	# Plot labels and legends
	plot!(legend=:topright, title="Absolute error as function of x", 
		xlabel="range (0,2) ", ylabel="Absolute Error of Central difference")

	printstyled("Saving Aboslute Error part (a) graph\n")
	savefig(p, "./output/Abs_Error_cdf.png")
	
	println("\nAnswering part b")
	printstyled("The largest error occured at: $(error[2]) at x = 1\n")

	println("\nAnswering part c")
	printstyled("The smallest error was: $(error[1]) and $(error[2]) at x = 0 and x = 2\n")
end

function sick()

	print("Solving part (a)\n")

	sick = Vector{Int}([25, 26, 28, 28, 33, 40, 52, 58, 65, 73, 74, 93, 105, 132,
                    144, 156, 164, 186, 210, 230, 239, 254, 268, 284, 317, 349,
                    408, 455, 489, 514, 568, 619, 675, 716, 780, 814, 829, 873,
                    914, 950])

	days = Vector{Int}(1:40)

	central = Float64[]

	second_central = Float64[]

	# Loops over 2-39 which is same as 1-38
	for i in 2:(length(sick)-1) 
		# Using central differenciation
		cdf = ( sick[i+1] - sick[i-1] ) / (days[i+1] - days[i-1])
		push!(central, cdf)
	end

	p = plot(days, sick, label="people sick", lc=:black, lw=2)
	
	plot!(days[2:end-1], central, label="First derivative", lc=:blue, lw=2)

	# Plot labels and legends
	plot!(legend=:topright, title="Number of population sick + 1st derivative", 
		xlabel="Number of days", ylabel="Number of people sick")

	printstyled("Saving Sick population and first derivative graph\n")
	savefig(p, "./output/Sick_1stDerivative.png")

	for i in 2:(length(sick)-1)
		# Using central differenciation for second derivative
		cdf_sec = ( sick[i+1] - 2*sick[i] + sick[i-1]) / ((days[i+1] - days[i])^2)
		push!(second_central, cdf_sec)
	end

	p = plot(days, sick, label="people sick", lc=:black, lw=2)
	plot!(days[2:end-1], second_central, label="Second derivative", lc=:blue, lw=2)

	# Plot labels and legends
	plot!(legend=:topright, title="Number of population sick + 2nd derivative", 
		xlabel="Number of days", ylabel="Number of people sick")

	printstyled("Saving Sick population and second derivative graph\n")
	savefig(p, "./output/Sick_2ndDerivative.png")

	println("\nAnswers to part c")
	printstyled("Sick people on day 14: $(sick[14])\n")
	printstyled("Sick people (first derivative) on day 14: $(central[14])\n")
	printstyled("Since the derivative at day 14 is positive the number of people getting sick will increase.")
	printstyled(" About 12 people will get sick per day since derivative is the slope. \n")
	printstyled("Assuming to have 396 people sick.\n")
	printstyled("Sick people on day 36: $(sick[36])\n")
	printstyled("The estiminate is way lower than the acutal number. It is not within the 25% of actual value.\n")

	println("\nSolving part (d)")

	log_sick = log.(sick)

	log_central = Float64[]

	log_second_der = Float64[]

	# Loops over 2-39 which is same as 1-38
	for i in 2:(length(sick)-1) 
		# Using central differenciation
		cdf = ( log_sick[i+1] - log_sick[i-1] ) / (days[i+1] - days[i-1])
		push!(log_central, cdf)
	end

	p = plot(days, log_sick, label="people sick", lc=:black, lw=2)
	
	plot!(days[2:end-1], log_central, label="First derivative", lc=:blue, lw=2)

	# Plot labels and legends
	plot!(legend=:topright, title="Number of population sick + 1st derivative", 
		xlabel="Number of days", ylabel="log of number of people sick")

	printstyled("Saving log of sick population and first derivative graph\n")
	savefig(p, "./output/Log_Sick_1stDerivative.png")

	for i in 2:(length(sick)-1)
		# Using central differenciation for second derivative
		cdf_sec = ( log_sick[i+1] - 2*log_sick[i] + log_sick[i-1]) / ((days[i+1] - days[i])^2)
		push!(log_second_der, cdf_sec)
	end

	p = plot(days, log_sick, label="people sick", lc=:black, lw=2)
	plot!(days[2:end-1], log_second_der, label="Second derivative", lc=:blue, lw=2)

	# Plot labels and legends
	plot!(legend=:topright, title="Number of population sick + 2nd derivative", 
		xlabel="Number of days", ylabel="log of number of people sick")

	printstyled("Saving log of sick population and second derivative graph\n")
	savefig(p, "./output/log_Sick_2ndDerivative.png")

	println("\nAnswers to part e")
	printstyled("Sick people on day 14: $(log_sick[14])\n")
	printstyled("Sick people (first derivative) on day 14: $(log_central[14])\n")
	printstyled("Since the derivative at day 14 is positive the number of people getting sick will increase.")
	printstyled(" About 0.08352 people will get sick per day since derivative is the slope.\n")
	printstyled("Assuming to have about 6.71 people sick.\n")
	printstyled("Sick people on day 36: $(log_sick[36])\n")
	printstyled("The estiminate very close to the real value. It is within the 25%\n")

end


function different_derivatives()

	# log(x) = ln(x) in julia
	function f(x)
		return sec(log(x))
	end

	function f_prime(x)
		return (1/x)*sec(log(x))*tan(log(x))
	end

	x = 0.5:0.1:4.5
	h = 0.1

	printstyled("Computing forward difference\n")
	forward_diff = Float64[]
	forward_error = Float64[]

	for i in x
		fwd = ( f(i+h) - f(i) ) / h
		push!(forward_diff, fwd)

		error = (fwd - f_prime(i)) / f_prime(i)
		push!(forward_error, error)
	end

	printstyled("Computing backward difference\n")

	backward_diff = Float64[]
	backward_error = Float64[]

	for i in x
		bwd = ( f(i) - f(i-h) ) / h 
		push!(backward_diff, bwd)
		error = (bwd - f_prime(i)) / f_prime(i)
		push!(backward_error, error)
	end

	printstyled("Computing central difference\n")

	central_diff = Float64[]
	central_error = Float64[]

	for i in x
		cdf = ( f(i+h) - f(i-h) ) / (2h)
		push!(central_diff, cdf)
		error = (cdf - f_prime(i)) / f_prime(i)
		push!(central_error, error)
	end

	p = plot(x, forward_error, label="Forward Diff error", lc=:black, lw=2)
	
	plot!(x, backward_error, label="Backward Diff error", lc=:blue, lw=2)
	
	plot!(x, central_error, label="Central Diff error", lc=:green, lw=2)

	# Plot labels and legends
	plot!(legend=:topright, title="Errors of different approximations", 
		xlabel="x-axis", ylabel="y-axis")

	printstyled("Saving graphs of errors of different approximations\n")
	savefig(p, "./output/approximations_errors.png")

end

function partial_derivatives()

	h = 0.25

	∂fₓ = (1.0625-1.8125)/(2h)

	printstyled("Answer to part (a): $(∂fₓ)\n")

	∂f_y = (0.5625 - 2*0.375 - 0.3125)/(h^2)

	printstyled("Answer to part (b): $(∂f_y)\n")

	∂f_xy = (1.3125-0.5625-1.8125+1.0625)/4*(h^2)

	printstyled("Answer to part (c): $(∂f_xy)\n")

end


# Same as if __name__ == '__main__'
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
