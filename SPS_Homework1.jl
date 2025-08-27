# SPS 5100 Research Techniques - Homework 1
using Plots

# Will put this in function when I have time


# Problem 1
########################################
printstyled("\nSolving Problem Number 1\n")

numbers = Int64[]

for i in range(length=10)
	push!(numbers, i)
	end

printstyled("List of Numbers: $(numbers)\n")
########################################

# Problem 2
########################################
printstyled("\nSolving Problem Number 2\n")

letters = [Char(c) for c in 'a':'a'+7]
printstyled("First 8 characters: $(letters)\n")

########################################

# Problem 3
########################################
printstyled("\nSolving Problem Number 3\n")

print("Enter an integer: ")
number = readline()

#Converting string to int
num = parse(Int, number)

if num % 2 == 0
	printstyled("The number $(num) is even\n")
else
	printstyled("The number $(num) is odd\n")
end

########################################

# Problem 4
########################################
printstyled("\nSolving Problem Number 4\n")

even_num = Int64[]
odd_num = Int64[]

for i in numbers
	if i % 2 == 0
		push!(even_num, i)
	else
		push!(odd_num, i)
	end
end

printstyled("Even numbers are: $(even_num)\n")
printstyled("Odd numbers are: $(odd_num)\n")

########################################

# Problem 5
########################################

printstyled("\nSolving Problem Number 4\n")

# Data points given in the assignment
x = [0.503, 1.283, 3.155, 4.353, 6.042]
y = [0.894, 0.284, -1.170, -0.373, 1.058]


# make a sin(x) plot from 0 to 2π 
θ = range(0, 2π, length=200)
f = sin.(θ)
p = plot(θ, f, label="sin(x)", lc=:black, lw=2)

# Plotting the dots above
scatter!(x, y, label="data", mc=:red)

# Plot labels and legends
plot!(legend=:bottomright, title="Sine with data points", 
	xlabel="x", ylabel="y")

printstyled("saving graph\n")
savefig(p, "sine.png")
########################################

# Problem 6
########################################

printstyled("\nSolving Problem Number 6\n")

# Part a) - Second half of the list
half = letters[1:4]
printstyled("Solution to Problem 6a: $(half)\n")

# Part b) - All entried whose index is odd
odd_letters = letters[1:2:end]
printstyled("Solution to Problem 6b: $(odd_letters)\n")

########################################

