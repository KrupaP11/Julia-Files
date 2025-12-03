# SPS 5100 Research Techniques - Homework 5

using Plots
gr()
using Statistics
using LinearAlgebra
include("Question4_HW5.jl")

function main()

	# Creating output folder
	# specify the directory name
	dir = "./output/HW5"

	# check if it exists, if not, create it
	if !isdir(dir)
	    mkpath(dir)
	    println("Directory created: $dir")
	end

	# Question 1
	println("Solving Question 1:")
	Question1()

	# Question 2
	println("\nSolving Question 2:")
	Question2()

	# Question 3
	println("\nSolving Question 3:")
	Question3()

	# Question 4
	println("\nSolving Question 4:")
	Question4()

	# Question 5
	println("\nSolving Question 5 - part a:")
	Question5_part_a()
	
	println("\nSolving Question 5 - part b:")
	Question5_part_b()

	# Question 6
	println("\nSolving Question 6:")
	Question6()

end

# Creating the functions before solving questions

# Newton's Method for question 2 and 4
function Newton(f, f′, x_start, max_iter=200, tol=1e-8)
	
	x_n = x_start

	for i in 1:max_iter

		# Getting the y-value
		f_x = f(x_n)

		# Checking if f'(x) is zero.
		if f′(x_n) == 0
			return nothing
		end

		# Checking if we are within the error bounds
		if abs(f_x) < tol
			return x_n, i
		end

		# If we are not within the limit we gotta updated our guess
		# Implementing Newton's method
		x_n -= f(x_n) / f′(x_n)

		#println(i, " ", x_n,)
	end

	print("Exceeded iteration limit without solution")
	return nothing
end

# Bisection Method for question 4
function Bisection(f, a, b, max_iter=200, tol=1e-8)

	# Compute f at the bottom of the interval
    f_a = f(a)

    # Print an initial update statement
    #println("0 ", (b+a)/2, " ", (b-a)/2)
    x_new = NaN

    # Main bisection method loop
    for i in 1:max_iter

        # Update interval midpoint and function value
        x_mid = (a+b) / 2.0
        f_mid = f(x_mid)

        # If f at midpoint has opposite sign to f at a, then it is the new b
        if f_a * f_mid < 0
            b = x_mid
        # Otherwise the midpoint replaces a
        else
            a, f_a = x_mid, f_mid
        end

        # Compute the midpoint of the updated interval and compare to
        #   the previous midpoint to see if we're close enough
        x_new = (a+b)/2.0
        xdiff = abs(x_new-x_mid)
        #println(i, " ", x_new, " ", xdiff)

        if abs(xdiff/x_new) < tol
        	return (x_new, i)
            break
        end
    end
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

function Question1()

	C = [
		0 1 0 0 0;
		0 0 1 0 0;
		0 0 0 1 0;
		0 0 0 0 1;
		243.433 -191.793 -23.954 30.319 0.832;
	]
	
	println("Part a) \nThe matrix C is: ")
	display(C)

	n = size(C, 1)
    Q_total = I(n)

    for k in 1:1000
        # Built-in QR decomposition
        Q, R = qr(C)
        Q = Matrix(Q)   # Convert to a dense matrix
        C = R * Q       # Form next C
        Q_total *= Q    # Accumulate eigenvector matrix

        # Convergence test: lower triangle → ~ 0
        # tril returns the lower-triangular part of
        if norm(tril(C, -1)) < 1e-10
            break
        end
    end

    # The diagonal entries of A ≈ eigenvalues
    λ = diag(C)
    println("Part b) \nAll eigenvalues aka roots are: $(λ)")
end

function Question2()

	# Writing the function and derivative of the function 
	
	function f(x)
		return (1+x^2)^(-1) - 1/2
	end

	function f′(x)
		return (-2*x)/(1+2*x^2+x^4)
	end

	# Range of x-values.
	x = -3:0.01:3

	# Plot to check what the function looks like
	p = plot([i for i in x], [f(i) for i in x], label="f(x)", lw=2, lc=:blue, xlabel="x values", ylabel="f(x) values")

	# Implementing Newton's Method

	# plot the roots over the older functions
	root = Float64[]
	x_pos = Float64[]
	x_neg = Float64[]

	for x_start in x
		xᶜ = Newton(f, f′, x_start)

		if xᶜ === nothing
			push!(root, x_start)
			elseif xᶜ == 1
				push!(x_pos, x_start)
			elseif xᶜ == -1
				push!(x_neg, x_start)
			else
				push!(root, x_start)
		end
	end

	plot!(legend=:bottomright)

	scatter!([root], [1], markershape=:cross, markeralpha=1, markersize=1, label="root ≢ ±1", markerstrokecolor=:red)

	scatter!([x_pos], [1], markershape=:circle, markersize=8, label="root = 1", markerstrokecolor=:blue)
		
	scatter!([x_neg], [1], markershape=:utriangle, markersize=8, label="root = -1", markerstrokecolor=:black)

	savefig(p, "./output/HW5/Newton's_Method_plot.png")
    println("Saved ./output/HW5/Newton's_Method_plot.png")
end

function Question3()

	#=
		For bisection method we take midpoint
		length of interval get cut in half every time. 
		l = (b-a)/2^n
		We want the interval to be less than 10^-16 so l = 10^-16

		10^-16 = (1-0)/2^n

		Now we solve for n

		10^16 = 2^n
	=#

	n = log2(10^16)
	println("The number of iterations to reach ϵ≈1e-16 is: $(n)")
end

function Question4()

	# Writing the function
	# The root according to desmos is x = 0.7032
	function f(x)
		return ℯ^x - 2 - 0.01/(x^2) + 2e-6/(x^3)
	end

	function f′(x)
		return ℯ^x + 2(0.01)/x^3 - 3(2e-6)/(x^4)
	end

	# Writing out the intervals
	I₁ = [0.01, 3]
	I₂ = [0.01, 9]
	I₃ = [0.01, 27]
	I₄ = [0.01, 81]

	# Need to show for each interval how many evaluations it takes to get the root

	# Implementing Bisection method
	I₁_x_bis, I₁_n_bis = Bisection(f, 0.01, 3, 200, 1e-6)
	I₂_x_bis, I₂_n_bis = Bisection(f, 0.01, 9, 200, 1e-6)
	I₃_x_bis, I₃_n_bis = Bisection(f, 0.01, 27, 200, 1e-6)
	I₄_x_bis, I₄_n_bis = Bisection(f, 0.01, 81, 200, 1e-6)


	println("\nRoot and iteration per interval for bisection method:")
	@info  "Interval = $(I₁)" I₁_x_bis I₁_n_bis
	@info  "Interval = $(I₂)" I₂_x_bis I₂_n_bis
	@info  "Interval = $(I₃)" I₃_x_bis I₃_n_bis
	@info  "Interval = $(I₄)" I₄_x_bis I₄_n_bis

	# Implementing Newton's method

	# will need midpoint per interval
	m₁ = (3+0.01)/2
	I₁_x_new, I₁_n_new = Newton(f, f′, m₁, 200, 1e-6)

	m₂ = (9+0.01)/2
	I₂_x_new, I₂_n_new = Newton(f, f′, m₂, 200, 1e-6)

	m₃ = (27+0.01)/2
	I₃_x_new, I₃_n_new = Newton(f, f′, m₃, 200, 1e-6)

	m₄ = (81+0.01)/2
	I₄_x_new, I₄_n_new = Newton(f, f′, m₃, 200, 1e-6)

	println("\nRoot and iteration per interval for Newton's method:")
	@info  "Interval = $(I₁)" I₁_x_new I₁_n_new
	@info  "Interval = $(I₂)" I₂_x_new I₂_n_new
	@info  "Interval = $(I₃)" I₃_x_new I₃_n_new
	@info  "Interval = $(I₄)" I₄_x_new I₄_n_new

	# Implementing Ridders' Method

	I₁_x_rid, I₁_n_rid = ridders(f, 0.01, 3, 200, 1e-6)
	I₂_x_rid, I₂_n_rid = ridders(f, 0.01, 9, 200, 1e-6)
	I₃_x_rid, I₃_n_rid = ridders(f, 0.01, 27, 200, 1e-6)
	I₄_x_rid, I₄_n_rid = ridders(f, 0.01, 81, 200, 1e-6)


	println("\nRoot and iteration per interval for Ridder's method:")
	@info  "Interval = $(I₁)" I₁_x_rid I₁_n_rid
	@info  "Interval = $(I₂)" I₂_x_rid I₂_n_rid
	@info  "Interval = $(I₃)" I₃_x_rid I₃_n_rid
	@info  "Interval = $(I₄)" I₄_x_rid I₄_n_rid
end

function Question5_part_a()

	println("\nWorking on part a)")

	#defining the constants
	c = 3*10^8
	h = 6.626*10^(-34)
	kB = 1.380649*10^(-23)

	function I(λ, T)

		# lambda.^(-5) dot is needed because we are treating lambda as vector.\
		β = (h*c) ./ (λ .* kB .* T)

		#Implement an if statement for e^50
		# Since \beta is a vector cause \lambda is a vector need to handle this differently.

		# For β > 50, exp(β) ≈ huge, so denominator ≈ exp(β)
		denom = @. ifelse(β > 50, exp(β), exp(β) - 1)
		return 2 * π * h * c^2 .* λ.^(-5) ./ denom
	end

	# Determining the range of temperature in Kelvin
	# adding units cause of the package
	T = 2700:100:10000

	# Using Wein's displacement law: \lambda_max = b/T
	#	We get the range for wavelengths

	#λ = 200:100:2000 # in linear space[nm]
	# Wavelength range (in meters)
    λ = range(1e-7, 3e-6, length=10000)

	# matrix for lambda max
	λ_max = Float64[]

	for i in T

		# Going through all the wavelengths
		Intensity = I.(λ, i)

		# Get the lambda max for that temperature
		λ_max_T = λ[argmax(Intensity)]

		push!(λ_max, λ_max_T)
	end

	T_log = log10.(T)
	λ_log = log10.(λ_max)

	# TODO Change the axis labels
	p = plot(T_log, λ_log, lw=2, lc=:blue)
    plot!(title="Temperature vs. Wavelength log plot", xlabel="log(Temperature) [K]", ylabel="log(λ_max) [m]")

    # Build the design matrix [x 1]
	X = [T_log ones(length(T_log))]
	
	# Solve least squares: β = [slope, intercept]
	β = X \ λ_log
	slope, intercept = β[1], β[2]
	λ_log_fit = intercept .+ slope .* T_log
	
	println("Slope = ", slope)
	println("Intercept = ", intercept)
	println("Estimated b = ", 10^intercept, " m·K")

	plot!(T_log, λ_log_fit, lw=1, lc=:red, ls=:dash, label="Fitted curve, slope=$(round(slope, digits=9))")

    savefig(p, "./output/HW5/Weins_Displacement.png")
    println("Saved ./output/HW5/Weins_Displacement.png\n")
end

function Question5_part_b()

	println("\nWorking on part b)")

	#defining the constants
	c = 3*10^8
	h = 6.626*10^(-34)
	kB = 1.380649*10^(-23)

    λ₁ = 390e-9
    λ₂ = 750e-9

    # Coding x^3/(e^x - 1)
    #f(x) = x^3 / (exp(x) - 1)

    function f(x)
    	# For x > 50, e^(x) ≈ huge, so denominator ≈ e^(β)
		denom = @. ifelse(x > 50, exp(x), exp(x) - 1)
		return x^3 ./ denom
    end

    function η(T)
    	x₁ =  h * c/(λ₂ * kB * T)
    	x₂ = h * c/(λ₁ * kB * T)

    	I = simpson(f, x₁, x₂, 10000)
    	return I * 15/(pi^4)
    end

    function golden_ratio(f, a, b, tol=1e-1)
    	ϕ = (sqrt(5) + 1) / 2
    	c = b - (b-a) / ϕ
    	d = a + (b-a) / ϕ

    	f_c = f(c)
    	f_d = f(d)

    	while abs(b - a) > tol
    		if f_c > f_d
    			a = c
    			c = d
    			f_c = f_d
    			d = a + (b - a)/ϕ
    			f_d = f(d)
    		else 
    			# f(c) < f(d)
    			b = d
    			d = c
    			f_d = f_c
    			c = b - (b - a)/ϕ
    			f_c = f(c)
    		end
    	end

    	if f_c > f_d

    		return (a+c)/2, f((a+c)/2)
    	else
    		return (b+d)/2, f((b+d)/2)
    	end

    end

	T_best, η_best = golden_ratio(η, 2700, 10000)

	println("Temperature for maximum efficiency: $(T_best)")
	println("Visible efficiency at the that temperature: $(η_best)")

	T = 2700:100:10000

	p = plot(T, [η(i) for i in T], lw=2, lc=:blue, xlabel="Temperature [K]", ylabel="Efficiency of Blackbody",
		title="Efficiency vs. Temperature")

	savefig(p, "./output/HW5/Efficiency_Temperature.png")
    println("Saved ./output/HW5/Efficiency_Temperature.png\n")

end

function Question6()
	
	# Implementing given function
	function f(x₀, x₁)
		return x₀^2 - 2*x₀ + x₁^4 - 2x₁^2 + x₁
	end

	∇f(x₀, x₁) = [2*x₀ - 2;
				4*x₁^3 - 4*x₁ + 1]

	println("\nWorking on part a)")
	println("The partial derivative with respect to x₀ is: 
		\n 2x₀-2 = 0 ⟹  x₀ = 1")
	println("The partial derivative with respect to x₁ is: 
		\n 4x₁^3 - 4x₁ + 1 ⟹  x₁ = {-1.1072, 0.2696, 0.8376}")

	println("Working on part b)")

	x₀_range = 0:0.5:2
	x₁_range = -1.5:1:1.5

	Z = [f(x₀, x₁) for x₁ in x₁_range, x₀ in x₀_range]

	#X0 = repeat(reshape(x₀, 1, :), length(x₁), 1)  
	#X1 = repeat(reshape(x₁, :, 1), 1, length(x₀))

	# Getting the output values from the function
	#Z = f.(X0, X1)

	p = plot(
        x₀_range, x₁_range, Z,
        st = :surface,
        c = :plasma,
        alpha = 0.8,         # make it slightly transparent
        line_z = Z,          # adds wireframe lines
        lw = 0.5,
        xlabel = "x₀",
        ylabel = "x₁",
        zlabel = "f(x₀, x₁)",
        title = "Surface Plot of f(x₀, x₁)",
        zlims=(-2, 3),
        camera = (60, 30)
    )

    savefig(p, "./output/HW5/2-D_Graph.png")
    println("Saved ./output/HW5/2-D_Graph.png")

    println("Working on part c)")

    I_1 = [1.0, 0.26]
    I_2 = [1.0, 0.27]

    # Gradient descent method

    function gradient_descent(a, b; γ=0.01, tol=1e-8, max_iter=1000)
    	x = [a, b]
    	path = [copy(x)]

    	for i in 1:max_iter
    		g = ∇f(x[1], x[2])
    		x_new = x .- γ .* g
    		push!(path, copy(x_new))

    		# if the new x_value and the old are within the error
    		#	or if the gradient is within the error return the path
        	if norm(x_new .- x) < tol || norm(g) < tol
            	return path
        	end
        	x = x_new
    	end
    	return path
	end

	path_a = gradient_descent(1.0, 0.26)
	path_b = gradient_descent(1.0, 0.27)

	#z_starts = [f(I_1[1], I_1[2]), f(I_2[1], I_2[2])]

	# Convert paths to arrays for plotting
	#xa = [p[1] for p in path_a]; ya = [p[2] for p in path_a]; za = f.(xa, ya)
	#xb = [p[1] for p in path_b]; yb = [p[2] for p in path_b]; zb = f.(xb, xb)

	xa, ya = [p[1] for p in path_a], [p[2] for p in path_a]
    xb, yb = [p[1] for p in path_b], [p[2] for p in path_b]
    za, zb = f.(xa, ya), f.(xb, yb)
	
	# Overlay paths (lines + points). Use different markers/colors.
	plot3d!(xa, ya, za; lw=2, marker=:circle, ms=3, label="path a (1.0, 0.26)", color=:red)
	plot3d!(xb, yb, zb; lw=2, marker=:diamond, ms=3, label="path b (1.0, 0.27)", color=:green)

	# also mark starting points and final points
	scatter3d!([1.0, 1.0], [0.26, 0.27], [f(1.0, 0.26), f(1.0, 0.27)];
        marker=:star5, ms=8, color=:yellow, label="starts")
    scatter3d!([last(path_a)[1], last(path_b)[1]],
        [last(path_a)[2], last(path_b)[2]],
        [f(last(path_a)[1], last(path_a)[2]), f(last(path_b)[1], last(path_b)[2])];
        marker=:x, ms=8, color=:blue, label="finals")

	savefig(p, "./output/HW5/2-D_Graph_path.png")
    println("Saved ./output/HW5/2-D_Graph_path.png")
end

# Same as if __name__ == '__main__'
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
