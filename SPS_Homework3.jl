# SPS 5100 Research Techniques - Homework 3

using Plots
using Printf

function main()
	# Creating output folder
	# specify the directory name
	dir = "./output/HW3"

	# check if it exists, if not, create it
	if !isdir(dir)
	    mkpath(dir)
	    println("Directory created: $dir")
	end

	# Question 1
	println("Solving Question 1:")
	#Simp_Trap()

	# Question 2
	println("\nSolving Question 2:")
	Mid_Trap()

	println("\nSolving Question 4:")
	Space_Object()

	println("\nSolving Question 5:")
	Romberg_Error_Analysis()

	println("\nSolving Question 3:")
	println("So the computer keeps killing the program. I don't know what I am doing wrong.")
	#Electric_potential()
end

# Question 1
function Simp_Trap()

	function f(x)
		return (sin(x))^4
	end

	a = 0
	b = 2π

	# Arrays for appending

	x_values = Vector{Vector{Float64}}()
	Trap_sum = Float64[]
	Simp_sum = Float64[]
	trap_rel_error = Float64[]
	simp_rel_error = Float64[]

	# This is the panel size. 
	# Ranging from 4 to 20 in order to get 5 to 21 points.
	n = 4:2:20

	for i in n
		
		# h is the panel width.
		h = (b-a)/i
		
		# Number of points or x-values
		# not i-1 because Julia starts counting from 1 not 0
		x = [a + j*h for j in 0:i]
		push!(x_values, x)

		# Getting the y-values
		y_values = [f(x) for x in x]

		# Calculating the Trapeziod Rule
		########################################################
		
		# Computing the first and the last points
		# "end" is for getting the last element of the array
		T_sum = h*( (1/2)*y_values[1] + sum(y_values[2:end-1]) + (1/2)*y_values[end] )
		push!(Trap_sum, T_sum)
	
		# Relative error: (computated - true)/true
		rel_error = abs((T_sum - 2.3561944901923449) / 2.3561944901923449)
		push!(trap_rel_error, rel_error)

		# Calculating the Simpson's Rule
		########################################################

		# Since n is already even here we don't have to worry. We can still check
		if i%2 != 0
			throw(ArgumentError("Simpson's Rule requires an even number of panels."))
		end

		# Getting all the coefficients first and then will multiply with f(x) values
		# vcat stands for vertical concatenation. So vcat([1, 2], [3, 4]) -> [1, 2, 3, 4]
		# the question mark is a form of if-else statement
			# condition ? value_if_true : value_if_false
		coefficients = vcat([1], [j%2 == 1 ? 4 : 2 for j in 2:i], [1])

		# In julia .* is like python's numpy * 
		S_sum = (h/3) * sum(coefficients .* y_values)
		push!(Simp_sum, S_sum)
		
		rel_error = abs((S_sum - 2.3561944901923449) / 2.3561944901923449)
		push!(simp_rel_error, rel_error)

	end

	# The Answer to question 1 part a
	println("Part a: \n")
	println("The answer using Trapezoid Rule is: $(Trap_sum[1])")
	println("The relative using Trapezoid Rule is: $(trap_rel_error[1])")
	println("The answer using Simpson's Rule is: $(Simp_sum[1])")
	println("The relative using Simpson's Rule is: $(simp_rel_error[1])")

	function parabola(x, x0, x1, x2, y0, y1, y2)
		ℓ0 = ((x - x1)*(x - x2))/((x0 - x1)*(x0 - x2))
    	ℓ1 = ((x - x0)*(x - x2))/((x1 - x0)*(x1 - x2))
    	ℓ2 = ((x - x0)*(x - x1))/((x2 - x0)*(x2 - x1))
    	return y0*ℓ0 + y1*ℓ1 + y2*ℓ2
	end

	# Plotting

	## I used ChatGPT for this. 
	## I spent more time fighting Julia for plot than getting the data for the plot -__-

	x5 = x_values[1]          
    y5 = f.(x5)

    xf = range(a, b, length=400)
    p = plot(xf, f.(xf), label="f(x)=sin^4(x)", lw=2)

    # nodes and connected lines (trapezoid edges)
    plot!(x5, y5, lw=2, marker=:circle, lc=:black, label="nodes (5 pts)")

    # draw trapezoid panels as filled shapes (semi transparent)
    for j in 1:length(x5)-1
        xpoly = [x5[j], x5[j+1], x5[j+1], x5[j], x5[j]]
        ypoly = [0.0, 0.0, y5[j+1], y5[j], 0.0]
        plot!(xpoly, ypoly, seriestype=:shape, fillalpha=0.15, linealpha=0.0, label=false)
    end

    # Simpson: two parabolas: (x0,x1,x2) and (x2,x3,x4)
    xs1 = range(x5[1], x5[3], length=200)
    ys1 = [parabola(xx, x5[1], x5[2], x5[3], y5[1], y5[2], y5[3]) for xx in xs1]
    plot!(xs1, ys1, lw=2, lc=:red, label="Simpson panel 1")

    xs2 = range(x5[3], x5[5], length=200)
    ys2 = [parabola(xx, x5[3], x5[4], x5[5], y5[3], y5[4], y5[5]) for xx in xs2]
    plot!(xs2, ys2, lw=2, lc=:red, label="Simpson panel 2")

    plot!(legend=:topright, title="Trapezoid vs Simpson panels (5 points)", xlabel="x", ylabel="y")
    savefig(p, "./output/HW3/Simp_Trap_Panel.png")
    println("Saved ./output/HW3/Simp_Trap_Panel.png")

    # Part (c) error plot
    points = collect(5:2:21)   # number of points
    p2 = plot(points, trap_rel_error, lw=2, marker=:circle, label="Trapezoid")
    plot!(points, simp_rel_error, lw=2, marker=:square, label="Simpson")
    plot!(legend=:topright, title="Relative error vs points", xlabel="Points", ylabel="Relative error (log10)")
    savefig(p2, "./output/HW3/Simp_Trap_Error.png")
    println("Saved ./output/HW3/Simp_Trap_Error.png")
	
end

# Question 2
function Mid_Trap()

	true_value = 0.6061244734187699
	
	a = 0
	b = 3

	function f(x)
		return sin(ℯ^x)
	end

	# For Simpson's Rule
	function simpson(f, a, b, n)
		
		# Checking if n is odd if add 1 to make it even
		if n % 2 == 1
			n += 1
		end

		# Computing the sum
		h = (b - a)/n
		x = range(a, stop=b, length=n+1)
		y = f.(x)
		coefficients = vcat([1], [i%2 == 1 ? 4 : 2 for i in 2:n], [1])
		
		return (h/3) * sum(coefficients .* y)
	end

	# For Midpoint Rule
	function midpoint(f, a, b, n)
		h = (b-a)/n
		x_mid = a .+ ((0:n-1) .+ 0.5) .* h
		return h*sum(f.(x_mid))
	end

	# Look for the number of points that matches the requirement
	function find_points(method; error=1e-5, max_n =10_000_000)
		n = 100
		while n<= max_n
			approx_val = method(f, a, b, n)
			rel_error = abs((approx_val - true_value) / true_value)

			if rel_error <= error
				return n, rel_error
			end
			
			n += 10
			#print(n)
		end
		return nothing, nothing
	end

	println("Part a: \n")
	mid_n, mid_error = find_points(midpoint)
	println("Midpoint rule needs $(mid_n+1) number of points.")
	println("So the answer is about 630 points")
	println("The relative error is: $(mid_error)")

	println("Part b: \n")
	simp_n, simp_error = find_points(simpson)
	println("Simpson's rule needs $(simp_n+1) number of points.")
	println("So the answer is about 29500 points")
	println("The relative error is: $(simp_error)")

end

# Question 3
function Electric_potential()
	
	# Define constants
	ϵₒ =  8.854187e-12 #[F/m]
	λ = 10e-6 #[C/m]

	a = -1
	b = 1
	rel_acc = 10e-6

	function f(s)
		r = sqrt((s^2) + (s^4 - 3)^2)
		return 1/r		
	end

	# I chose simpson's rule because it is easier to code
	function simpson(f, a, b, n)
		
		# if odd add one
		if n%2 == 1
			n+=1
		end

		h = (b - a)/n

		x = range(a, stop=b, length=n+1)
		y = f.(x)
		
		coefficients = vcat([1], [i%2 == 1 ? 4 : 2 for i in 2:n], [1])
		
		return (h/3) * sum(coefficients .* y)
	end

	# using richardson error estimation

	function richardson_simp(f, a, b, rel_acc)
		n = 100
		Iₙ = simpson(f, a, b, n)
		while true
			n *= 2
			I_2n = simpson(f, a, b, n)
			error = abs(I_2n - Iₙ)/15
			relative_error = error/abs(I_2n)

			if relative_error < rel_acc
				return I_2n, relative_error, n
			end
			I_n = I_2n
			print(n)
		end
	end

	# Integral value
	I_val, relative_error, n_used = richardson_simp(f, a, b, rel_acc)
	V = λ/(4*π*ϵₒ) * I_val

	println("Relative Error Estimation is: $(relative_error)")
	println("Number of panels: $(n_used)")
	println("Electric_potential V = $(V)")
end

# Question 3

function Space_Object()
	
	# Constants
	H₀ = 67.66 #[km/s/Mpc]
	c = 3e5 #[km/s]
	Ωᵣ = 0.4165/(H₀^2)
	Ωₖ = 0
	Ωₘ = 0.3111 - 0.5*Ωᵣ
	Ωᵥ = 0.6889 - 0.5*Ωᵣ

	# z goes from 0 to 10
	a = 0
	b = 10
	n = 200

	z_range = range(a, stop=b, length=200)

	function f(z)
		E = sqrt(Ωᵣ*(1+z)^4 + Ωₘ*(1+z)^3 + Ωₖ*(1+z)^2 + Ωᵥ)
		return 1/E
	end

	# Trapezoid rule
	function trap(f, a, b, n)

		h = (b - a)/n
		x = range(a, stop=b, length=n+1)
		y = f.(x)
		return h*( (1/2)*y[1] + sum(y[2:end-1]) + (1/2)*y[end] )
	end

	# Richardson extrapolation
	function Richardson(Iₙ, I_2n)
		return I_2n + ( (I_2n - Iₙ)/3 )
		# p = 2 for trapezoid rule. p^2 - 1
	end

	trap_Dc = zeros(length(z_range))
	trap_2n_Dc = zeros(length(z_range))
	Rich_Dc = zeros(length(z_range))
	
	for (i, z) in enumerate(z_range)
	    if z == 0
	        trap_Dc[i] = 0.0
	        trap_2n_Dc[i] = 0.0
	        Rich_Dc[i] = 0.0
	    else
	        In = trap(f, 0, z, n)
	        I2n = trap(f, 0, z, 2*n)
	        Ir = Richardson(In, I2n)
	
	        trap_Dc[i] = (c/H₀)*In
	        trap_2n_Dc[i] = (c/H₀)*I2n
	        Rich_Dc[i] = (c/H₀)*Ir
	    end
	end


	p = plot(z_range, Rich_Dc, label="Richardson extrapolated", lw=2, lc=:blue)
	plot!(z_range, trap_2n_Dc, label="Trapezoid (2n panels)", lw=2, lc=:red)
    plot!(legend=:bottomright, title="Comoving distance vs. Redshift", xlabel="Redshift z", ylabel="Comoving distance D_C [Mpc]")

    savefig(p, "./output/HW3/comoving_redshift.png")
    println("Saved ./output/HW3/comoving_redshift.png")

end

function Romberg_Error_Analysis()

	a = 0
	b = π

	true_value = 0.64400474375804016

	ns = [1+2^i for i in 2:12]

	# Romberg Table
	R = zeros(length(ns), length(ns))

	# Arrays for errors
	romberg_errors = Float64[]			# part (a)
	expected_errors = Float64[]			# part (b)
    richardson_errors = Float64[]       # part (c)

    h = (b-a)/(ns[1]-1)
    x = [a + i*h for i in 0:(ns[1]-1)]

    function f(x)
    	return sin(ℯ^x)
    end

    # Function for trap rule

    function trap_rule(f, points, h)
    	total = sum(f.(points[2:end-1]))
    	total += 0.5*(f.(points[1]) + f(points[end]))
    	return h*total
    end

    # Adaptive Trap rule
    function trap_rule_adaptive(f, points, h)
    	total = sum(f.(points[2:2:end-1]))
    	return h*total
    end

    trap_adap = trap_rule(f, x, h)
    R[1, 1] = trap_adap
    push!(romberg_errors, abs(R[1, 1]- true_value))

    # Building the romberg table

    for (k, n) in enumerate(ns[2:end])
    	h *= 0.5
    	x = [a + i*h for i in 0:n]
    	trap_adap *+ 0.5
    	trap_adap += trap_rule_adaptive(f, x, h)

    	i = k+1
    	R[i, 1] = trap_adap

    	for j in 2:i
            R[i,j] = (4^j * R[i,j-1] - R[i-1,j-1]) / (4^j - 1)
    	end

    	# Part (a): actual Romberg error (bottom-right of row)
        push!(romberg_errors, abs(R[i,i] - true_value))

        # Part (b): expected error using R[i-1,j] and R[i-1,j-1]
        # Formula: (R[i-1,j] - R[i-1,j-1]) / (4^j - 1)
        j = i  # current diagonal
        expected_err = j > 2 ? abs((R[i-1,j-1] - R[i-1,j-2]) / (4^(j-1) - 1)) : NaN
        push!(expected_errors, expected_err)
		
		# Part (c): single-step Richardson extrapolation error using R[i,0] and R[i-1,0]
        richardson_val = (4*R[i,1] - R[i-1,1]) / 3.0
        richardson_err = abs(richardson_val - true_value)
        push!(richardson_errors, richardson_err)
    end

    # Plotting

    p = plot(ns, romberg_errors, marker=:circle, label="Actual romber error", lc=:blue)
    plot!(ns[2:end], expected_errors, marker=:square, label="Expected error", lc=:red)
    plot!(ns[2:end], richardson_errors, marker=:cross, label="Richardson Step error", lc=:black)

    plot!(xlabel="Number of points used", ylabel="Error (absolute)", title="Romberg Integration error analysis", legend=:bottomright)

    savefig(p, "./output/HW3/Romberg_error_plot.png")
    println("Saved ./output/HW3/Romberg_error_plot.png")

end



# Same as if __name__ == '__main__'
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
