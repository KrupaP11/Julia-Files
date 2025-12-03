# SPS 5100 Research Techniques - Homework 6

using Plots

function main()

	# Creating output folder
	# specify the directory name
	dir = "./output/HW6"

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

end


function Question2()

	# Coding in the two time-dependent populations
	#	alpha is the per-capita growth rate of prey
	#	beta is the rate at which interaction reduce the prey population
	#	gamma is the per-caputa death rate of pedator
	#	delta is the rate at which interaction increase the predator population

	# The predator
	function model(t)
		h, f = t
		dhdt = α*h - β*f*h     	#Predator
		dfdt = -γ*f + δ*f*h     #Prey
		return [dhdt, dfdt]	
	end

	println("\nWorking on part a)")

	# Initial population of foxes
	f0 = 4000

	# Initial population of hares
	h0 = 35000

	# Parameters for prey
	α = 0.2
	β = 2.5E-5
	# Parameters for predator
	γ = 0.15
	δ = 1.5E-5

	# Deciding time
	# T = 2*pi / sqrt(gamma*alpha) = 36.3 years 

	# Two cycles would be 
	t_final = 80.0  # years 
	dt = 0.01 		# Time step
	N = Int(floor(t_final/dt)) + 1
	t = range(0, stop=t_final, length=N)

	# RK4 method
	y = zeros(N, 2)
	y[1, :] = [h0, f0]

	for i in 1:(N-1)
		yi = @view y[i, :]
    	k1 = model(yi)
    	k2 = model(yi .+ 0.5*dt*k1)
    	k3 = model(yi .+ 0.5*dt*k2)
    	k4 = model(yi .+ dt*k3)
    	y[i+1, :] = yi .+ dt * (k1 .+ 2k2 .+ 2k3 .+ k4) / 6
    	y[i+1, :] = max.(y[i+1, :], 0.0)   # avoid negative values
	end

	hares = y[:, 1]
	foxes = y[:, 2]

	p = plot(t, hares, label="Hares (prey)", xlabel="Time (years)", ylabel="Population")
	plot!(t, foxes, label="Foxes (predator)", title="Lotka–Volterra Model (RK4)")

	savefig(p, "./output/HW6/Lotka-Volterra_Model.png")
    println("Saved ./output/HW6/Lotka-Volterra_Model.png")

    println("\nWorking on part b)")

    # Find local maxima of hares
	peaks_idx = []
	for i in 2:(N-1)
	    if hares[i] > hares[i-1] && hares[i] > hares[i+1]
	        push!(peaks_idx, i)
	    end
	end

    # Cycle length (time between first two maxima)
	if length(peaks_idx) >= 2
	    cycle_length = t[peaks_idx[2]] - t[peaks_idx[1]]
	    println("Cycle length ≈ $(cycle_length) years")
	else
	    println("Not enough peaks to determine cycle length.")
	end

    println("\nWorking on part c)")

    # Global maximum hare population
	imax = argmax(hares)
	max_hare = hares[imax]
	max_time = t[imax]
	println("Maximum hare population: $(round(max_hare, digits=3))")
	println("hares at t ≈ $(round(max_time, digits=3)) years")
end

function Question3()
	
	# The four model equations

	function model(S, E, I, R)
	
		dSdt = (-β*I*S)/N
		dEdt = (β*I*S)/N - σ*E
		dIdt = σ*E - γ*I
		dRdt = γ*I

		return [dSdt, dEdt, dIdt, dRdt]
	end

	#Paramters
	# Transimission efficiency
	β = 0.514
	# rate at which exposed individuals become infectious
	σ = 0.50
	# the clearance rate; recovery rate
	γ = 0.20

	println("\nWorking on part a)")
	
	# Deciding time
	t_final = 300
	dt = 1
	t = 0:1:300
	n_steps = length(t)

	# Preallocate arrays
	S = zeros(Float64, n_steps)
	E = zeros(Float64, n_steps)
	I = zeros(Float64, n_steps)
	R = zeros(Float64, n_steps)

	# susceptible individuals
	S0 = 6E6 - 2

	# infectious individuals
	I0 = 2

	R0 = 0.0

	E0 = 0.0

	S[1], E[1], I[1], R[1] = S0, E0, I0, R0
	# Total population
	N = S0 + E0 + I0 + R0

	# RK4 method

	for n in 1:(n_steps-1)
		s, e, i, r = S[n], E[n], I[n], R[n]
    	k1s, k1e, k1i, k1r = model(s, e, i, r)
    	k2s, k2e, k2i, k2r = model(s + 0.5*dt*k1s, e + 0.5*dt*k1e, i + 0.5*dt*k1i, r + 0.5*dt*k1r)
    	k3s, k3e, k3i, k3r = model(s + 0.5*dt*k2s, e + 0.5*dt*k2e, i + 0.5*dt*k2i, r + 0.5*dt*k2r)
    	k4s, k4e, k4i, k4r = model(s + dt*k3s,     e + dt*k3e,     i + dt*k3i,     r + dt*k3r)
	
    	S[n+1] = s + dt*(k1s + 2.0*k2s + 2.0*k3s + k4s)/6.0
    	E[n+1] = e + dt*(k1e + 2.0*k2e + 2.0*k3e + k4e)/6.0
    	I[n+1] = i + dt*(k1i + 2.0*k2i + 2.0*k3i + k4i)/6.0
    	R[n+1] = r + dt*(k1r + 2.0*k2r + 2.0*k3r + k4r)/6.0

    	#push!(S, S[n+1])
    	#push!(E, E[n+1])
    	#push!(I, I[n+1])
    	#push!(R, R[n+1])
	end

	p1 = plot(t, S, label="Susceptible", xlabel="Days", ylabel="Population",
          title="SEIR populations over 0–300 days", lw=2, lc=:red)
	plot!(t, R, label="Recovered", lc=:blue, lw=2)
	plot!(t, I, label="Infected", lc=:black, lw=2)
	plot!(t, E, label="Exposed", lc=:green, lw=2)

	savefig(p1, "./output/HW6/SEIR_population.png")
    println("Saved ./output/HW6/SEIR_population.png")

    println("\nWorking on part b)")

    peak_idx = argmax(I)            # index in Julia (1-based)
	peak_day = t[peak_idx]
	peak_value = I[peak_idx]
	
	println("Peak infectious day (integer): ", Int(round(peak_day)))
	println("Peak infectious value: ", round(peak_value; digits=1), " individuals")

	println("\nWorking on part c)")
	println("Attached as a pdf")

	println("\nWorking on part d)")

	# Using κ(t) = ln[I(t)/I(t-1)] for t = 1:300 (discrete consecutive days)
	# align indices with t; leave kappa[1] = NaN for t=0
	κ = fill(NaN, n_steps)  

	# Start at 2 because we need 0, and 1 to start 
	#	if we start at 1 then there isn't a value before 0
	for j in 2:n_steps
	    if I[j-1] > 0.0 && I[j] > 0.0
	        κ[j] = log(I[j] / I[j-1]) / (t[j] - t[j-1])   # denominator is 1.0 here
	    else
	        κ[j] = NaN
	    end
	end

	p2 = plot(t[2:end], κ[2:end], marker=:circle, markersize=3, lw=1,
          xlabel="Days (t)", ylabel="κ (growth rate per day)",
          title="Estimated daily growth rate κ for days 1:300")

	savefig(p2, "./output/HW6/growth_rate.png")
    println("Saved ./output/HW6/growth_rate.png")

end

function Question4()

	println("\nPart a is attached as pdf.")
	
	println("\nWorking on part b)")

	global funcevals = 0


	# Coding the ODE
	function f(x, y)
		global funcevals
		funcevals +=1
		dydx = 1/2 * (x-1)*(2-y)
		return dydx
	end

	# Establishing state vector
	function fvec(t, S)
		global funcevals
		funcevals += 1
	
		# Unpack the state vector
		x = S[0]
		y = S[1]
		vx = S[2]
		vy = S[3]
	
		# Compute the derivatives
		fx = vx
		fy = vy
		fvx = -2*vx*V
		fvy = -2*vy*V
	
		return Float64[fx, fy, fvx, fvy]
	end

	# Implementation of Runge-Kutta-Fehlberg method for solving 1D 1st-
	#	order ODEs.

	function rkf45(f, a, b, y0, tol=1e-6)
		# Very crude initial guess for h
	
		h = (b-a)/100.0
		x = a
		y = y0
	
		while true
			# Adssume we don't exit this step; change if needed
			break_now = false
	
			# Check whether current step size would carry past b,
			#	and adjust h if so
			if x+h > b
				h = b - x
				break_now  = true
			end
	
			# Perform the RK step
			k0 = h*f(x,y)
			k1 = h*f(x+h/4, y+k0/4)
			k2 = h*f(x + 3*h/8, y + 3*k0/32 + 9*k1/32)
			k3 = h*f(x + 12*h/13, y + 1932*k0/2197 - 7200*k1/2197 + 7296*k2/2197)
			k4 = h*f(x+h, y + 439*k0/216 - 8*k1 + 3680*k2/513 - 845/4104*k3)
			k5 = h*f(x+h/2, y - 8*k0/27 + 2*k1 - 3544*k2/2565 +
				1859*k3/4104 - 11*k4/40)
			ystar = y + 25*k0/216 + 1408*k2/2565 + 2197*k3/4104 - k4/5
			ynew = y + 16*k0/135 + 6656*k2/12825 + 
			28561*k3/56430 - 9*k4/50 + 2*k5/55
	
			# Check the error; update x, y if allowed
			err = abs(ystar - ynew)
			if err<tol
				x += h
				y = ystar
			end
	
			# Exit if instucted to
			if break_now
				break
			end
	
			# Adjust step size
			if err != 0.0
				h *= 0.8*(tol/err)^0.2
			else
				h *= 0.5
			end
	
		end
		return y
	end

	# Implementation of 4th-order Runge-Kutta method for solving
	#	simultaneous 1st-order ODEs.
	
	function rk4_vec(f, a, b, N, S0)
		h = (b - a)/Float64(N)
		ts = a .+ (0:N) .* h
	
		# Make sure this has the right shape to hold the state vectors
		Ss = zeros(Float64, (N, 4))
		S = S0
	
		for (i, t) in enumerate(ts)
			Ss[i] = S
			k0 = h*f(t, S)
			k1 = h*f(t+0.5*h, S+0.5*k0)
			k2 = h*f(t+0.5*h, S+0.5*k1)
			k3 = h*f(t+h, S+k2)
			S += (k0 + 2*k1 + 2*k2 + k3)/6.0
		end
		return ts, Ss
	end
	
	# True solution
	ans(x) = 2 - ℯ^((x/2) - x^2/4)

	# Set initial quantities
	y0 = 1

	# Starting point of interval
	a = 0

	# Ending point of interval
	b = 5

	# Create the arrays to hold y-values

	ys_rk45 = []
	fs_rk45 = []
	rk45_errs = []

	# Fill out the array of tolerances to use for RKF45
	t = range(-1, -13, length=50)
	tolerance = 10 .^ t

	for tol in tolerance
		funcevals = 0
		push!(ys_rk45, rkf45(f, a, b, y0, tol))
		push!(fs_rk45, funcevals)
		push!(rk45_errs,  abs(ys_rk45[end] - ans(b)) / abs(ys_rk45[end]) )
	end

	p = scatter(fs_rk45, rk45_errs, xlabel="Function Evals", ylabel="Error", xscale=:log10, yscale=:log10, label="Errors")
	scatter!(fs_rk45, tolerance, label="Tolerance", xscale=:log10, yscale=:log10)

	savefig(p, "./output/HW6/rkf45_error.png")
    println("Saved ./output/HW6/rkf45_error.png")
end

# Global Constants
const T_lowest = 1E-8
const ρ_lowest = 1E-12
const ρ_min = 1E-6 

function Question5()
		
	println("Unfortunately, a lot of help from ChatGPT was needed in order to finish this problem.")

	# Solar units
	Msun = 1.98847e30       #[kg]
	Rsun = 6.96e8           #[m]
	Lsun = 3.83e26          #[W]

	# Physical Constants
	c = 3.0E8 							 #[m/s]
	G = 6.67430E-11         			 #[m^3 kg^-1 s^-2]
	k_B = 1.380649E-23				     #[m^2 kg s^-2 K^-1]
	α = 7.565E-16                        #[J m^-3 K^-4]

	# Star parameters
	Mstar = 100.0 * Msun    # total mass (100 M_sun)

	# Return mass density rho from Pressure and temperature
	# The restrictions ensures finite positive output
	function density_EOS(P::Float64, T::Float64)

		# Key problem we can't let P < 1/3 * a * T^4 because then density is negative
		P_rad = 1/3 * α * T^4

		ρ = (9.91E-28) * (P - P_rad) / (k_B * max(T, T_lowest))

		if ρ <= 0.0 || !isfinite(ρ)
			return ρ_min
		end

		return max(ρ, ρ_min)
	end

	# The restrictions allow a large finite positive value is the formula misbehaves.
	function opacity(ρ::Float64, T::Float64)

		ρ_safe = max(ρ, ρ_lowest)
		T_safe = max(T, T_lowest)

		# Making sure that the Temperature doesn't turn negative
		κ = 0.035 + 6.44E18 * (ρ/T_safe^(3.5))         #[m^2 kg^-1]

		# Keeping kappa finite and positive

		if !isfinite(κ) || κ <= 0 
			return 1E6
		end
		return κ
	end

	function epsilon(ρ::Float64, T::Float64)

		ρ_safe = max(ρ, ρ_lowest)
		T_safe = max(T, T_lowest)

		T6 = T_safe*1E-6  #[K]

		T₆ = max(T6, 1E-12)

		# avoid throwing DomainError by using safe values
		ϵ = 0.136 * ρ * T₆^(-2/3) * ℯ^(-33.80*T₆^(-1/3)) + 
		4.69E11 * ρ^2 * T₆^-3 * ℯ^(-4403/T₆)                 #[W kg^-1]

		if !isfinite(ϵ) || ϵ < 0
			return 0.0
		end

		return ϵ
	end

	# Coding the ODEs
	# Given enclosed mass Mᵣ and state vecotr y = [r, P, Lr, T] return 
	# 		the derivatives with respect to enclosed mass Mᵣ
	function derivatives(Mᵣ::Float64, y::Vector{Float64})

		r, P, Lᵣ, T = y

		ρ = density_EOS(P, T)

		# if radius is very small or density is invalid we will set the derivatives
		#	and avoid calling opcaity and epsilon

		if r <= 0.0 || !isfinite(r) || ρ <= 0.0 || !isfinite(ρ)
	        return Float64[0.0, 0.0, 0.0, 0.0]
	    end

    	# Now it is safe to call opacity and epsilon
		κ = opacity(ρ, T)
		ϵ = epsilon(ρ, T)

		dLᵣdMᵣ = ϵ
		drdMᵣ = 1.0 / (4.0 * π * r^2 * ρ)
		dPdMᵣ = - (G * Mᵣ) / (4.0 * π * r^4)

		# restriction against zero/negative temperature
		if T <= 0.0 || !isfinite(T)
			dTdMᵣ = 0.0
		else
			denom = 64.0 * π^2 * α * c * T^3 * r^4
			if denom == 0.0 || !isfinite(denom)
				dTdMᵣ = 0.0
			else
				dTdMᵣ = - (3 * κ * Lᵣ) / denom
			end
		end

		#gotta make sure that nothing is NaN or infinite after the computation

		if any(x -> (!isfinite(x)), (drdMᵣ, dPdMᵣ, dLᵣdMᵣ, dTdMᵣ))
			return Float64[0.0, 0.0, 0.0, 0.0]
		end

		return Float64[drdMᵣ, dPdMᵣ, dLᵣdMᵣ, dTdMᵣ]
	end

	println("\nWorking on part a)")

	function rk4(Mᵣ::Float64, y::Vector{Float64}, h::Float64)
		k1 = derivatives(Mᵣ, y)
		k2 = derivatives(Mᵣ + 0.5 * h, y + 0.5*h .* k1)
    	k3 = derivatives(Mᵣ + 0.5 * h, y + 0.5*h .* k2)
    	k4 = derivatives(Mᵣ + h, y + h .* k3)
	    return y + (h/6.0) .* (k1 .+ 2 .* k2 + 2 .* k3 .+ k4)	
	end

	# Mass grid generator
	# Return an array of N+1 mass points from 0 to Mstar with concetration
	#	near the center controlled by exponent p (>1 concetrates centrally)

	function make_mass_grid(Mstar::Float64, N::Int; p::Float64=3.0)

		#undef tells Julia to allocate the memory for the array 
		#	but to skip the step of initializing its elements.
		Mgrid = Vector{Float64}(undef, N+1)
		for i in 0:N
        	frac = (i / N)^p
        	Mgrid[i+1] = frac * Mstar
    	end
    	return Mgrid
	end

	# Integrator
	"""
	Integrate outward in enclosed mass from the center to Mstar using a 
	non-uniform mass grid. The very first non-zero mass shell is initialized
	consistently using Mᵣ = Mgrid[2] and the cnetral density rhoc

	Returns Mᵣ_arr, Ys_arr, positive_flag
	
	Mᵣ_arr : vector of enclosed mass points (from first non-zero shell to Mstar)
	Ys_arr: matrix with rows [r, P, Lr, T] at each Mᵣ_arr point
	positive_flag: true if integration did not produce nonpositive pressure or NaN/Inf
	"""

	function integrate_star_debug(Pc::Float64, Tc::Float64, N::Int; debug_steps::Int=8, pgrid::Float64=3.0)
		
		# Check the radiation pressure at the center
		P_rad_c = (1.0/3.0) * α * Tc^4
		if Pc <= P_rad_c
        	#println("SKIP: Pc=$(Pc) <= P_rad(Tc)=$(P_rad_c). density_EOS will return ρ_min.")
        	r0 = 1E-6
        	y0 = Float64[r0, Pc, 0.0, Tc]
        	return [0.0], reduce(vcat, (reshape(y0, 1, :) for _ in 1:1)), false
    	end

    	# Create mass grid
    	Mgrid = make_mass_grid(Mstar, N; p=pgrid)
		
		# Initialize at firs non-zero mass point
		if length(Mgrid) < 2
			error("Mass grid must have at least 2points")
		end

		ρ_c = density_EOS(Pc, Tc)

		M1 = Mgrid[2]
		# Set the initial radius by inverting M = (4/3) π r^3 ρ  => r = (3 M / (4 π ρ))^(1/3)
    	r1 = (3.0 * M1 / (4.0 * π * ρ_c))^(1.0/3.0)

    	L1 = epsilon(ρ_c, Tc) * M1

    	# initial state at M = M1
    	y = Float64[r1, Pc, L1, Tc]

    	# storage
    	Mᵣ_val = Float64[]
    	Ylist = Vector{Vector{Float64}}()

    	positive = true

    	# intergrating from index 2
    	for i in 2:(length(Mgrid)-1)
    		Mᵣ = Mgrid[i]
    		h = Mgrid[i+1] - Mgrid[i]

    		# debug prints for first few steps
        	if i-1 <= debug_steps
        	    P_rad = (1.0/3.0) * α * y[4]^4
        	    ρ_dbg = density_EOS(y[2], y[4])
        	    κ_dbg = opacity(ρ_dbg, y[4])
        	    ϵ_dbg = epsilon(ρ_dbg, y[4])
        	    drdMᵣ, dPdMᵣ, dLdMᵣ, dTdMᵣ = derivatives(Mᵣ, y)
        	    #println("step=$(i-1)  M_r=$(Mᵣ)  r=$(y[1])  P=$(y[2])  T=$(y[4])")
        	    #println("       P_rad=$(P_rad)  rho=$(ρ_dbg)  kappa=$(κ_dbg)  epsilon=$(ϵ_dbg)")
        	    #println("       dr/dMᵣ=$(drdMᵣ)  dP/dMᵣ=$(dPdMᵣ)  dL/dMᵣ=$(dLdMᵣ)  dT/dMᵣ=$(dTdMᵣ)")
        	end
        	
        	# Not letting anything go off to infinity
        	if !(isfinite(y[1]) && isfinite(y[2]) && isfinite(y[4])) || y[2] <= 0.0 || y[1] <= 0.0 || y[4] <= 0.0
            	positive = false
            	#println("Pre-step sanity failed at Mᵣ=$(Mᵣ). y=$(y)")
            	break
        	end

        	# Take RK4 step in mass coordinate
       	 	y_new = rk4(Mᵣ, y, h)

       	 	# append
       	 	push!(Mᵣ_val, Mᵣ + h)
       	 	push!(Ylist, copy(y_new))

       	 	if y_new[2] <= 0.0 || any(x -> !isfinite(x), y_new)
       	 		positive = false
            	#println("Stopped: nonpositive P or NaN/Inf after stepping from Mr=$(Mr). P=$(y_new[2])")
       	 		break
       	 	end

       	 	# accept the step
       	 	y = y_new
       	 end

       	 if isempty(Mᵣ_val)
       	 	Mᵣ_arr = [M1]
       	 	Ys_arr = reduce(vcat, (reshape(Float64[r1, Pc, L1, Tc], 1, : ) for _ in 1:1))
       	 	return Mᵣ_arr, Ys_arr, false
       	 else
       	 	Mᵣ_arr = collect(Mᵣ_val)
       	 	Ys_arr = reduce(vcat, (reshape(v, 1, :) for v in Ylist))
       	 	return Mᵣ_arr, Ys_arr, positive
       	 end
	end

	# Debugging here
	function test_point(Pc::Float64, Tc::Float64, N::Int=200)
		println("\n=== Testing Pc=$(Pc), Tc=$(Tc) ===")
		M_arr, Ys_arr, ok = integrate_star_debug(Pc, Tc, N; debug_steps=6)
    	if !ok
    	    println("Integration finished: ok=$(ok). Final surface P = $(Ys_arr[end, 2])")
    	else
    	    println("Integration finished successfully. Final surface P = $(Ys_arr[end,2])")
    	end
	end

	function find_core_parameters(N)
	
	    # Sun's core temp is 10E6 K
	    T_grid = 10 .^ range(log10(5E6), log10(3E8), length=12)
	
	    # Sun's core pressure is 1.2E16 Pa
	    P_grid = 10 .^ range(15, 19, length=12)
	
	    best = nothing
	    best_metric = 1E99
	
	    fail_count = 0
	    explode_count = 0
	    total_tested = 0
	
	    for Tc in T_grid
	        for Pc in P_grid

				total_tested += 1
				Mᵣ_arr, Ys_arr, positive = integrate_star_debug(Pc, Tc, N)
	
	            if !positive
	                fail_count += 1
	                P_end = Ys_arr[end, 2]
	                #println("Failed at Pc=$(Pc), Tc=$(Tc).  Surface P = $P_end")
	                continue
	            end
	
	            # This is a valid integration
	            P_final = Ys_arr[end, 2]
	
	            if !isfinite(P_final)
	                explode_count += 1
	                println("Integration blew up: Pc=$Pc Tc=$Tc  P_final=$P_final")
	                continue
	            end
	
	            if P_final < best_metric
	                best_metric = P_final
	                best = (Pc, Tc, P_final, Mᵣ_arr, Ys_arr)
	            end
	        end
	    end
	
	    #println("Grid points tested: $(length(T_grid) * length(P_grid))")
	    #println("Failed integration count = $fail_count")
	    #println("Diverged (NaN/Inf) count = $explode_count")
	
	    if best === nothing
	        println("No valid coarse-grid solution found.")
	        return nothing
	    end
	
	    Pc0, Tc0 = best[1], best[2]
	    #println("\nBest coarse-grid point:")
	    #println("Pc0 = $Pc0   Tc0 = $Tc0   P_surface = $(best[3])")
	
	    # --- Refined local search ---
	    Pvals = range(0.5*Pc0, 1.5*Pc0, length=9)
	    Tvals = range(0.5*Tc0, 1.5*Tc0, length=9)
	
	    best2 = nothing
	    best_metric2 = 1E99
		
	    for Tc in Tvals
	        for Pc in Pvals
	            Mᵣ_arr, Ys_arr, positive = integrate_star_debug(Pc, Tc, N)
	
	            if !positive
	                continue
	            end
	
	            P_final = Ys_arr[end, 2]
	            if !isfinite(P_final)
	            	continue
	            end

	            if P_final < best_metric2
	                best_metric2 = P_final
	                best2 = (Pc, Tc, P_final, Mᵣ_arr, Ys_arr)
	            end
	        end
	    end
	
	    if best2 === nothing
	        #println("Fine grid failed — only coarse solution available.")
	        return best, best[4]  # coarse-grid result only
	    end
	
	    #println("Fine grid solution found!")
	    return best2, best2[4]
	end


	N = 100

	result_data = find_core_parameters(N)

	if result_data === nothing
	    println("Failed to find a positive surface pressure solution in the search grid.")
	else
		result, Mᵣ_arr = result_data
		Pc_best, Tc_best, P_final, Mᵣ_arr, Ys_arr = result
		r_arr = Ys_arr[:, 1]
		P_arr = Ys_arr[:, 2]
		L_arr = Ys_arr[:, 3]
		T_arr = Ys_arr[:, 4]
	
		R_star = r_arr[end]
		L_star = L_arr[end]	

		println("Central pressure Pc = $(round(Pc_best, sigdigits=3)) Pa")
		println("Central temperature Tc = $(round(Tc_best, sigdigits=3)) K")
		println("Surface pressure (final) = $(round(P_arr[end], sigdigits=3)) Pa")
		
	
		println("\nWorking on part b)")
		println("Star radius R = $(round(R_star, sigdigits=3)) [m] = $(round(R_star/Rsun, digits=3)) [R_sun]")
		println("Star luminosity L = $(round(L_star, sigdigits=3)) [W] = $(round(L_star/Lsun, sigdigits=3)) [L_sun]")
		
		println("\nTwo-significant-figure central values:")
		println("Pc ≈ $(round(Pc_best, sigdigits=3)) Pa")
		println("Tc ≈ $(round(Tc_best, sigdigits=3)) K")
	end

	# Quick manual tests (left commented so user can enable as needed)
	#test_point(1.0e15, 3.0e8)   # expected to SKIP due to P_rad >> Pc
	#test_point(1.0e19, 3.0e8)   # should attempt integration; use debug prints


	println("\nWorking on part c)")

	Mᵣ_plot = Mᵣ_arr ./ Msun

	p1 = plot(Mᵣ_plot, r_arr,
     xlabel="M_r (M_sun)",
     ylabel="r(M_r) [m] (log scale)",
     yscale=:log10,
     title="Radius vs Enclosed Mass",
     grid=:both,
     legend=false,
     lw=2)

	p2 = plot(Mᵣ_plot, T_arr,
     xlabel="M_r (M_sun)",
     ylabel="T(M_r) [K] (log scale)",
     yscale=:log10,
     title="Temperature vs Enclosed Mass",
     grid=:both,
     legend=false,
     lw=2)

	p3 = plot(Mᵣ_plot, L_arr,
     xlabel="M_r (M_sun)",
     ylabel="L(M_r) [W] (log scale)",
     yscale=:log10,
     title="Luminosity vs Enclosed Mass",
     grid=:both,
     legend=false,
     lw=2)

	savefig(p1, "./output/HW6/Radius_mass.png")
    println("Saved ./output/HW6/Radius_mass.png")

	savefig(p2, "./output/HW6/Temperature_mass.png")
    println("Saved ./output/HW6/Temperature_mass.png")

	savefig(p3, "./output/HW6/Luminosity_mass.png")
    println("Saved ./output/HW6/Luminosity_mass.png")

    println("\nI am not sure if any of these answers were right or not. This question was a hard one.")
end

# Same as if __name__ == '__main__'
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end