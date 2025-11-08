using LinearAlgebra, Random

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
function array_create(n, val)
    # Initialize the array and vector to be empty
    A = Array{Float64}(undef, n, n)
    v = Array{Float64}(undef, n)
    
    # Create an array of eigenvalues
    v[:] = collect(val:val+n-1)
    v .*= 0.5
    # println("Input eigenvals:", v)
    
    # Create the diagonal matrix of eigenvalues
    D = I(n) .* v
    
    # Create matrix of eigenvectors
    A = reshape(collect(val:val+n*n-1), n, n)
    A = sqrt.(A)
    
    """ As coded, qr_decomp uses the standard Gram-Schmidt algorithm, which
    introduces significant orthogonalization errors (see function
    test_qrdec). The numpy implementation uses the modified Gram-
    Schmidt algorithm, which is far more accurate.
    Make sure one of the following two lines is commented out."""
    
    # Q, R = qr_decomp(A)
    Q, R = qr(A)  # Julia’s built-in QR decomposition
    
    # Multiply out to get an array with the specified eigenvalues/vectors
    A = (transpose(Q) * D) * Q
    
    # Create a random vector for the constant vector in Ax = b
    v[:] = [rand() for i in 1:n]
    
    # And send things back to the calling function
    return A, v
end

"""
Function to perform the QR decomposition of a square array.
It is assumed that the input array is square. Pass non-square arrays
at your peril!

Input argument:
(1) A_in: array to be factored

Outputs:
(1) Q: orthogonal square matrix
(2) R: upper triangular square matrix
"""
function qr_decomp(A_in)
    N = size(A_in, 1)
    Aprime = copy(A_in)
    Q = zeros(N, N)
    R = zeros(N, N)
    
    # Outer loop over columns
    for j in 1:N
        # Inner loop for orthogonalization
        for i in 1:j-1
            R[i, j] = dot(Q[:, i], A_in[:, j])
            Aprime[:, j] .-= R[i, j] .* Q[:, i]
        end
        R[j, j] = norm(Aprime[:, j])
        Q[:, j] = Aprime[:, j] ./ R[j, j]
    end
    
    return Q, R
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
function qr_method(A_in, max_itrs::Int=100)
    A = copy(A_in)
    
    for i in 1:max_itrs
        Q, R = qr_decomp(A)
        A = R * Q
        # println(i, diag(A))
        
        #TODO: this would be a good place to check the off-diagonal 
        # elements of A
    end
    
    qreigvals = diag(A)
    return qreigvals
end

"""
By putting our main function inside this if statement, we can safely
import the module from other scripts without having this code execute
every time
"""
if abspath(PROGRAM_FILE) == @__FILE__
    #TODO: Create an array
    #TODO: Compute its eigenvalues with the QR method
    
    np_eigvals, np_eigvecs = eigen(A)
    
    #TODO: Compare against numpy’s built-in functions, which
    # were just calculated
    #TODO: Plot the results
end
