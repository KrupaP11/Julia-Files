# SPS_Homework 5 Question 4
# Dr. Warren gave this code in python. Converting it to julia
#	in order to make it work with my code.

# implement this function below


#=

f ridders(f,a,b,max_itrs=200,tol=1e-8):
# Ensure that x0 < x1
if a < b:
x0 = a
x1 = b
else:
x0 = b
x1 = a
# Compute the function at the two endpoints and make sure they bracket
# a root
f0 = f(x0)
f1 = f(x1)
if f0 == 0.0:
return x0
elif f1 == 0.0:
return x1
if f0*f1 > 0.0:
raise Exception("Starting points for Ridders’ rule must bracket root")
# Main Ridders’ method loop
for i in range(max_itrs):
# Find the midpoint of the interval
x2 = (x0+x1)/2.0
f2 = f(x2)
# Apply Ridders’ formula to get x3
x3 = x2 + (x1-x2)*(f2/f0)/sqrt((f2/f0)**2-(f1/f0))
f3 = f(x3)
print(i,x3)
# Check to see if we can return: we took a very small step or
# our latest guess evaluates to 0
if abs(x3-x2) < tol or f3 == 0.0:
return x3
# Bracket the root as tightly as possible
if f3*f2 < 0.0:
# The root lies between x2 and x3, so set that as the interval for
# the next iteration
x0 = x2; f0 = f2
x1 = x3; f1 = f3
else:
# f3 and f2 have same sign, so x3 and x2 don’t bracket the root.
# Check whether f0 or f1 is the other endpoint
if f3*f1 < 0.0:
# x3 and x1 bracket the root, so replace x0 with x3
x0 = x3; f0 = f3
else:
# x3 and x0 bracket the root, so replace x1 with x3
x1 = x3; f1 = f3
print("Exceeded iteration limit without solution")
return None
=#