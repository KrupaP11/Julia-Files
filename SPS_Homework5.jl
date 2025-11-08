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
	Question1()

	# Question 2
	println("\nSolving Question 2:")
	#Question2()

	# Question 3
	println("\nSolving Question 3:")
	#Question3()

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

function Question1()

	C = [
		0 1 0 0 0;
		0 0 1 0 0;
		0 0 0 1 0;
		0 0 0 0 1;
		-243.433 -191.793 -23.954 -30.319 -0.832;
	]
	
	println("Part a) \n The matrix C is: $(display(C))")
end

# Same as if __name__ == '__main__'
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
