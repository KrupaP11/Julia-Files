# SPS 5100 Research Techniques - Homework 4

using Plots
using Random
using LinearAlgebra

function main()

	# Creating output folder
	# specify the directory name
	dir = "./output/HW4"

	# check if it exists, if not, create it
	if !isdir(dir)
	    mkpath(dir)
	    println("Directory created: $dir")
	end

	# Question 1
	println("Question 1 is attached as pdf.")

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
	println("\nSolving Question 5:")
	Question5()

	# Question 6
	println("\nSolving Question 6:")
	Question6()

end

# Creating a generalized Forward Substitution function
#= Apply forward substitution to solve a lower triangular matrix,
      using the method provided in Newman, pg 219.
    Input arguments:
      (1) L: square array, already in lower triangular form
      (2) b: vector representing right side of equation Lx=b
    Output:
      (1) x: solution vector to set of equations
=#
function forward_substitution(A, b)
	
	N = length(b)
    x = zeros(eltype(b), N)

    for m in 1:N
		# sum over previously computed entries (sum of empty slice is 0)
        x[m] = (b[m] - sum(A[m, 1:m-1] .* x[1:m-1])) / A[m,m]
    end
    return x
end

# Creating a generalized Backward Substitution function
#=Apply backsubstitution to solve a lower triangular matrix,
      using the method provided in Newman, pg 219.
    Input arguments:
      (1) A: square array, already in lower triangular form
      (2) b: vector representing right side of equation Ax=b
    Output:
      (1) x: solution vector to set of equations
=#
function backward_substitution(A, b)
	
	N = length(b)
	x = zeros(eltype(b), N)

	for m in N:-1:1
        # MUST use m+1:N for upper-triangular part
    	x[m] = (b[m] - sum(A[m, m+1:N] .* x[m+1:N])) / A[m,m]
    end
	return x
end

# Creating a generalized LU_Decomp function
#= Function to compute LU decomposition of an array
    Input argument:
      (1) A: square array (must also be non-singular, but that's a
            discussion for a later day)
    Outputs:
      (1) L_out: lower triangular array with diagonal entries = 1
      (2) U_out: upper triangular array
=#
function LU_Decomp(A)

	# Compute the LU decomposition
	#=
	A = A^-1 = I
	let A^-1 be B
	Break A down into LU
	=#

	U_out = Float64.(A)
	N = size(U_out, 1)
	L_out = Matrix{Float64}(LinearAlgebra.I, N, N)  # identity matrix

	for j in 1:(N-1)
		for i in j+1:N
			coeff = U_out[i,j]/U_out[j,j]
			# update row i
			U_out[i, j:N] .-= coeff .* U_out[j, j:N]
			L_out[i, j] = coeff
		end
	end

	return L_out, U_out
end

#=
Function to create an array using square roots of sequential
     integers.  Used in Gezerlis textbook rather than "array_create"
     as defined elsewhere in this script.
   Also creates a vector by rescaling the 0th row of the array
   Input arguments:
     (1) n: size of the array to create
     (2) val: starting value for the integer range
   Outputs:
     (1) A: square array of sequential square roots
     (2) bs: vector of rescaled values from 0th row of A
=#

function array_create_gezerlis(n, val)
	
	# n^2 - 1 is to ensure the last element isn't counted.
	# collect does the same thing as np.arrange.

	A = reshape(collect(val:(val + n^2 - 1)), n, n)
   	A = sqrt.(A)
   	
   	bs = (A[1, :]).^(2.1)

   	return A, bs
end

#= Function to perform the QR decomposition of a square array.  It is
      assumed that the input array is square.  Pass non-square arrays
      at your peril!
    Input argument:
      (1) A_in: array to be factored
    Outputs:
      (1) Q: orthogonal square matrix
      (2) R: upper triangular square matrix
=#

function QR_decomp(A)

	N = size(A, 1)
	Aprime = copy(A)

	Q = zeros(Float64, N, N)
	R = zeros(Float64, N, N)

	for j in 1:N
		for i in 1:j-1
			#.= and .* keep everything in-place and broadcasted, 
			#	avoiding temporary arrays.
			R[i, j] = dot(Q[:, i], A[:, j])
			Aprime[:, j] .-= R[i, j] .* Q[:,i]
		
		end

		R[j, j] = norm(Aprime[:, j])
		Q[:, j] = Aprime[:, j] ./ R[j, j]
	end
	
	return Q, R
end

function Question2()
	
	# Matrix Given
	A = [1 1 -1 1; 1 2 1 1; 2 3 -1 3; 3 4 1 6]

	# Finding the inverse
	println("Getting the inverse and determinant of the matrix")

	function inverse_matrix(A)
		n = size(A, 1)
		L, U = LU_Decomp(A)

		inv_A = zeros(Float64, size(A))

		for i in 1:n
			e = zeros(Float64, n)
			e[i] = 1.0

			y = forward_substitution(L, e)
			x = backward_substitution(U, y)

			inv_A[:, i] = x
		end
		return inv_A
	end

	inv_A = inverse_matrix(A)
	println()
	println("Inverse of Matrix A: \n $(inv_A[:,1]) \n $(inv_A[:,2]) \n $(inv_A[:,3]) \n $(inv_A[:,4])")

	function determinant(A)
		L, U = LU_Decomp(A)

		# Since L is a lower triangular matrix the determinant of L is 1

		# finding the determinant of U

		det = 1.0

		for i in 1:size(U, 1)
			det *= U[i, i]
		end

		return det
	end

	det = determinant(A)
	println()
	println("Determinant of Matrix A = $(det)")
end

function Question3()
	 #=
		Kirchoff's law states: The algebraic sum of all 
			currents entering and exiting a node (or junction) is zero. 
	 =#

	A = [
    9.5  -2.5   0.0  -2.0   0.0;
   -2.5  11.0  -3.5   0.0  -5.0;
    0.0  -3.5  15.5   0.0  -4.0;
   -2.0   0.0   0.0   7.0  -3.0;
    0.0  -5.0  -4.0  -3.0  12.0
	]

	b = [12.0; -16.0; 14.0; 10.0; -30.0]

	L, U = LU_Decomp(A)

	x̃ = forward_substitution(L, b)
	x = backward_substitution(U, x̃)

	println("The currents in the circut are: $(x) Ampere")
	println("The current passing through the 5Ω resistor is: $(x[2]) Ampere")
end

function Question4()
	
	# So we have (A - \lamda * I) v = 0 
	# First get the A - \lamda * I matrix
	# Then multiply that matrix with vector v = 0.
	# Using LU Decomp figure out what V is

	A = [
	1 2 0 5 0;
	-2 4 -2 -2 2;
	0 -3 5 0 0;
	5 0 2 -1 -5;
	2 4 0 0 -1
	]

	A = Float64.(A)

	# Let C = A - \lamda * I
	# Lamda given
	λ = 3

	# Computing the matrix
	C = A - λ * LinearAlgebra.I

	# Now we have C*v = 0
	# Or we have C*v = \lamda*v

	# Creating a matrix with 0
	b = zeros(size(A, 1))

	function gauss_elim(A, b)
		n = size(b, 1)

		for j in 1:n-1

			for i in (j+1):n
				coeff = A[i,j] / A[j,j]
				A[i,j:end] .-= coeff .* A[j,j:end]
				b[i] -= coeff * b[j]
			end
		end	
		x = backward_substitution(A, b)
		return x
	end

	x = gauss_elim(C, b)
	println("Solving Part a: \nI tried it with LU decomposition and turns out it couldn't do it.")
	println("Tried it with Guassian elimination and it also failed. Lastly I tried it with the in-build function.")
	println("The eigenvector associated with λ = 3 using gauss_elim: \n $(x) \n")

	# Using inbuild 

	# Get the nullspace → eigenvectors corresponding to λ
	# nullspace(C) returns an orthonormal basis for the nullspace of C
	#	v such that Cv = 0

	v = nullspace(C)
	
	println("Eigenvector associated with λ = 3 using in-build function:")
	display(v)

	# Couldn't use Guassian elimination how we have it because would need to pivot
	# Using nullspace was easier and safer in terms of calculations.

	n = size(A, 1)
    Q_total = I(n)

    for k in 1:1000
        # Built-in QR decomposition
        Q, R = qr(A)
        Q = Matrix(Q)   # Convert to a dense matrix
        A = R * Q       # Form next A
        Q_total *= Q    # Accumulate eigenvector matrix

        # Convergence test: lower triangle → ~ 0
        # tril returns the lower-triangular part of
        if norm(tril(A, -1)) < 1e-10
            break
        end
    end

    # The diagonal entries of A ≈ eigenvalues
    λ = diag(A)
    println("All four eigenvalues are: $(λ)")

end

function Question5()

	ϵ = Float64[]
	ϵ_jl = Float64[]

	n = 2:100

	for i in n

		# Create a T array using the function
		# Let the starting value be 2
		T, bs = array_create_gezerlis(i, 2)

		Q, R = QR_decomp(T)

		error = norm(transpose(Q)*Q - LinearAlgebra.I)
		push!(ϵ, error)

		Q_jl, R = qr(T)
		error = norm(transpose(Q_jl)*Q_jl - LinearAlgebra.I)
		push!(ϵ_jl, error)


	end

	p = plot(collect(n), ϵ, label="Error in code", lw=2, lc=:blue)
	plot!(collect(n), ϵ_jl, label="Error in Julia", lw=2, lc=:red)
    plot!(legend=:topleft, title="Error in Orthogonalization of Q", xlabel="Iterations", ylabel="Error")

    savefig(p, "./output/HW4/Error_in_Q.png")
    println("Saved ./output/HW4/Error_in_Q.png")

end

# Functions needed for Question6
"""
Generate a square array with specified eigenvectors and eigenvalues,
and also a vector of uniform random numbers between 0 and 1.

Input arguments:
(1) n: integer size of array and vector
(2) val: integer used to set element values of the array and vector

Outputs:
(1) A: square n by n matrix
(2) b: vector of length n
"""
function array_create(n::Int, val::Int)
	# Initialize empty arrays
    A = Matrix{Float64}(undef, n, n)
    v = Vector{Float64}(undef, n)

    # Create an array of eigenvalues
    v[:] = collect(val:val+n-1)
    v .*= 0.5

    # Create the diagonal matrix of eigenvalues
    D = I(n) .* v

    # Create matrix of eigenvectors
    A = reshape(collect(val:val+n*n-1), n, n)
    A = sqrt.(A)
	#= 
	As coded, qr_decomp uses the standard Gram-Schmidt algorithm, which
    introduces significant orthogonalization errors (see function
    test_qrdec). The numpy implementation uses the modified Gram-
    Schmidt algorithm, which is far more accurate.
    Make sure one of the following two lines is commented out.
    =#
    
    # Q, R = qr_decomp(A)
    # qr() is julia's built-in QR decomposition
    Q, R = qr(A)

    # Multiply out to get an array with the specified eigenvalues/vectors
    A = (transpose(Q) *D) *Q

    # Create a random vector for the constant vector in Ax = b
    v[:] = [rand() for i in 1:n]

    return A, v
end

"""
Function to compute the eigenvalues of a square array, via the
QR method.

There are no checks to ensure the array is square, so make sure
you never pass a non-square array or the results will be
unexpected.

Input arguments:
(1) A_in: the square array
(2) max_itrs (optional): the maximum number of times to
iterate through the QR method

Output:
(1) qreigvals: a vector holding the eigenvalues of A_in
"""

function qr_method(A, max_itrs)
	A = copy(A)

	for i in 1:max_itrs
		Q, R = QR_decomp(A)
		A = R * Q

		A_offdiag = abs.(A .- diag(A) .* I(size(A, 1)))
	end

	qreigvals = diag(A)
	return qreigvals
end

function Question6()

	A = Float64.(rand(1:10, 6, 6))

	qreigvals = qr_method(A, 100)

	eig_res = eigen(A)
	jl_eigvals = eig_res.values
	jl_eigvecs = eig_res.vectors

	# Compare with built-in eigenvalues
	# sorting to make sure they match
	error = abs.(sort(abs.(qreigvals)) .- sort(abs.(jl_eigvals)))
	println("Absolute errors between QR method and built-in eigenvalues: \n $(error)\n")
	
	p = bar(1:length(error), error,
	    label="Eigenvalue Error",
	    color=:orange,
	    xlabel="Index",
	    ylabel="Absolute Error",
	    title="Error between QR Method and Built-in Eigenvalues",
	    legend=:topright
	)

    savefig(p, "./output/HW4/QR_Method_error.png")
    println("Saved ./output/HW4/QR_Method_error.png")

end

# Same as if __name__ == '__main__'
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
