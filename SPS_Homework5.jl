# SPS 5100 Research Techniques - Homework 5

using Plots
using LinearAlgebra

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
	println("Solving Question 3:")
	#Question1()

	# Question 2
	println("\nSolving Question 2:")
	#Question2()

	# Question 3
	println("\nSolving Question 3:")
	Question3()

	# Question 4
	println("\nSolving Question 4:")
	#Question4()

	# Question 5
	println("\nSolving Question 5:")
	#Question5()

	# Question 6
	println("\nSolving Question 6:")
	#Question6()

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
			return x_n
		end

		# If we are not within the limit we gotta updated our guess
		# Implimenting Newton's method
		x_n -= f(x_n) / f′(x_n)
	end

	print("Exceeded iteration limit without solution")
	return nothing
end

# Bisection Method for question 4
function Bisection(f, a, b, max_iters = 200, tol=1e-8)

	# Compute f at the bottom of the interval
    f_a = f(a)

    # Print an initial update statement
    println("0 ", (b+a)/2, " ", (b-a)/2)
    x_new = NaN

    # Main bisection method loop
    for i in 1:max_itrs

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
        println(i, " ", x_new, " ", xdiff)

        if abs(xdiff/x_new) < tol
            break
        end
    end

    return x_new
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

	# Implimenting Newton's Method

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
	println("The number of iterations to reach ≈1e-16 is: $(n)")
end

function Question4()
	# Implement ridder's function first.

	#= need to call bisection, newton's and ridder's method and 
		show how many function evaluations does it take to get error
		within 1e-6.
	=#
end


# Same as if __name__ == '__main__'
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
