# SPS 5100 Research Techniques - Homework 1
using Plots
using BenchmarkTools

function main()

	# Creating output folder
	# specify the directory name
	dir = "./output"

	# check if it exists, if not, create it
	if !isdir(dir)
	    mkpath(dir)
	    println("Directory created: $dir")
	end

	# Question 1
	println("Solving question 1")
	Taylor_expansion_a(0.5)
	Taylor_expansion_b(10)

	# Question 2
	println("Solving question 2")
	Benchmarking()

	# Question 3
	println("Solving question 3")
	catastrophic_cancelation()

	# Question 4
	println("Solving question 4")
	checking_catastrophic_cancelation()

	# Question 5
	println("Solving question 5")
	plot_roundoff(100)
	

end


# Question 1 part a
function Taylor_expansion_a(x)

	# Putting an empty array of float64 to 
	#	collect the numbers
	sums = Float64[]

	# Collecting the errors
	ε_N = Float64[]

	# Defining terms
	value = 0.0

	# Putting in value
	exp_value = 1.6487212707001281
	
	# N iterations ranging from 0 to 20
	for i in 0:20
		
		# coding the equation 
		equ = x^i/factorial(i)

		# Adding to previous terms
		value += equ

		# Appending the values for plotting
		push!(sums, value)

		# Calculating the difference between 
		#	the expected value and calculated value
		error = abs(exp_value - value)

		# Appending the values for plotting
		push!(ε_N, log(error))

	end	

	# Plotting the absolute error as function of N

	# N-values
	N = range(0, 20)

	#log_ε = Float64[log(_) for _ in ε_N]

	p = plot(N, ε_N, label="Absolute error", lc=:black, lw=2)
	
	# Plot labels and legends
	plot!(legend=:topright, title="Absolute error as function of N", 
		xlabel="Number of Iterations", ylabel="log Absolute Error")

	printstyled("Saving Aboslute Error part (a) graph\n")
	savefig(p, "./output/Abs_Error_e^0.5.png")

end

# Question 1 part b
function Taylor_expansion_b(x)

	# Putting an empty array of float64 to 
	#	collect the numbers
	sums = Float64[]

	# Collecting the errors
	ε_N = Float64[]

	# Defining terms
	value = 0.0

	# Putting in value
	exp_value = 22026.465794806717
	
	# N iterations ranging from 0 to 20
	for i in 0:100
		
		# coding the equation 
		equ = x^i/factorial(big(i))

		# Adding to previous terms
		value += equ

		# Appending the values for plotting
		push!(sums, value)

		# Calculating the difference between 
		#	the expected value and calculated value
		error = abs(exp_value - value)

		# Appending the values for plotting
		push!(ε_N, log(error))

	end	

	# Plotting the absolute error as function of N

	# N-values
	N = range(0, 100)

	#log_ε = Float64[log(_) for _ in ε_N]

	p = plot(N, ε_N, label="Absolute error", lc=:black, lw=2)
	
	# Plot labels and legends
	plot!(legend=:topright, title="Absolute error as function of N", 
		xlabel="Number of Iterations", ylabel="log Absolute Error")

	printstyled("Saving Aboslute Error part (b) graph\n")
	savefig(p, "./output/Abs_Error_e^10.png")

end

# Question 2
function Benchmarking()

	# Size of the array
	Size = [50, 150, 500]

	# Empty array to hold time

	# This is for the brute force coded loop
	time_brute = Float64[]
	scale_brute = Float64[]
	# This is for the in-built matrix multiplication
	time_inbuilt = Float64[]
	scale_inbuilt = Float64[]

	for N in Size
		# Initializing the arrays
		A = rand(0:N, N, N)
		B = rand(0:N, N, N)
		C = zeros(N, N)

		# Timing the brute force multiplication
		t_loop = @belapsed begin
			# $ signs are needed because when using belapsed
			#	it doesn't recognize variables outside of it.
			#	The dollar sign helps but it as global variable
			for i in 1:$N, j in 1:$N, k in 1:$N
				$C[i, j] += $A[i, k]*$B[k, j]
			end # for i, j, k
		end # for time
			
		# Putting the time into array
		push!(time_brute, t_loop)

		# Reset C
		C .= 0
		# The dot is a broadcasting assignment
		# Doesn't assign a new memory allocation

		# Timing the in-built multiplication
		t_mul = @belapsed begin
			$C = $A * $B
		end # for time

		# Putting the time into array
		push!(time_inbuilt, t_mul)

	end # for N

	scale_inbuilt = log(time_inbuilt[3]/time_inbuilt[1]) / log(Size[3]/Size[1])

	scale_brute = log(time_brute[3]/time_brute[1]) / log(Size[3]/Size[1])

	## Testing ##
	printstyled("Time In-Built: $(time_inbuilt)\n")
	printstyled("Time Brute force: $(time_brute)\n")
	
	### TODO: Ask Dr. Warren if this is right or if 
	###			I need to get the slope of the start time and end time
	printstyled("Brute scaling: $(scale_brute)\n")
	printstyled("In built scaling: $(scale_inbuilt)\n")

	# Plotting the data above
	p = scatter(Size, time_inbuilt, label="In-Built Multiplication points", mc=:red)
	plot!(Size, time_inbuilt, label="In-Built Multiplication", lc=:red)
	scatter!(Size, time_brute, label="Brute Multication points", mc=:blue)
	plot!(Size, time_brute, label="Brute Multication", lc=:blue)

	# Plot labels and legends
	plot!(legend=:topleft, title="Time Comparison", 
		xlabel="Array size (N)", ylabel="Time [s]")

	printstyled("Saving Time Comparison graph\n")
	savefig(p, "./output/TimeComparison.png")
	
end

# Question 3
function catastrophic_cancelation()
	
	# Initializing varaibles
	x̃ = 12345678912345678
	ỹ = 12345678912345677

	# Calculating with int
	s = x̃^2 - ỹ^2
	printstyled("Using Ints: $(s)\n")

	# Rewriting as float
	x̃ = Float64(12345678912345678.0)
	ỹ = Float64(12345678912345677.0)

	# Calculating with float
	solution = x̃^2 - ỹ^2
	printstyled("Using Floats: $(solution)\n")

	# Using algebra to check
	s̃ = (x̃+ỹ)*(x̃-ỹ)
	printstyled("solution using multiplication: $(s̃) \n")

	#=
	The answer is closer to the float because both Python and Julia use
	arbitary precision to store integers, however, float requires a 
	8-bit precision which is causes loss of significant digits during 
	subtraction of nearly equal floats leading to catastrophic cancelation.
	=#

end

# Question 4
function checking_catastrophic_cancelation()
	x = 12345678912345678
	y = 12345678912345677

	# Part (a)
	ans_a = sqrt(x+1) - sqrt(x)
	printstyled("Answer without changes, part a: $(ans_a)\n")

	new_ans_a = 1/( sqrt(x+1) + sqrt(x) )
	printstyled("Answer after changes, part a: $(new_ans_a)\n")

	# Part (b)
	ans_b = 1/(x+1) - 2/x + 1/(x-1)
	printstyled("Answer without changes, part b: $(ans_b)\n")

	new_ans_b = 2/(x*(x^2-1))
	printstyled("Answer after changes, part b: $(new_ans_b)\n")

	# Part (c)
	ans_c = 1/sqrt(x) - 1/(sqrt(x-1))
	printstyled("Answer without changes, part c: $(ans_c)\n")

	new_ans_c = -1/( sqrt(x)*(x-1) + x*(sqrt(x-1)) )
	printstyled("Answer after changes, part c: $(new_ans_c)\n")

end

function roundoff(x, n)

	# This function was given to us.

	for i in 0:(n-1)
		x = sqrt(x)
	end

	for i in 0:(n-1)
		x = x^2
	end

	return x

end

function plot_roundoff(n)

	printstyled("Output of given code: \n")
	
	for x in (0.5, 5.0)
		printstyled( "$(x), $(roundoff(x, 100)) \n" )
	end

	# Getting the values
	y = [roundoff(5.0, i) for i in 0:n]
	#println(y)
	
	# Plotting	
	p = plot(0:n, y, label="f(5.0, $(n))", lc=:black, lw=2)
	
	# Plot labels and legends
	plot!(legend=:topright, title="Roundoff error", 
		xlabel="Number of Iterations", ylabel="Effect of roundoff errors on f(5.0, $(n))")

	printstyled("Saving Roundoff Error plot \n")
	savefig(p, "./output/roundoff_error.png")

end

# Same as if __name__ == '__main__'
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
