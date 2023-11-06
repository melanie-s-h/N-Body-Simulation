#========================================================================================#
"""
	NBodies
In this module a system of N interacting bodies moving in dim dimensions over time T in res timesteps
will be simulated.
N-body equations of motion are defined by: 
	F̲ᵢⱼ = mᵢ×(d²x̲ᵢ/dt²) = ∑ Gmᵢmⱼ×(x̲ⱼ-x̲ᵢ)×(1/|(x̲ⱼ-x̲ᵢ)|³), where the right side is summed up over N with j=1 & j≠i. 
These equations describe the motions of N mass points mᵢ, 
moving under the influence of their mutual attracting force given by Newton's law of gravitation.
with p̲ᵢ = mᵢ×(dx̲ᵢ/dt) follows a set of first order differential equations:
	ṗ̲ᵢ = dp̲ᵢ/dt = ∑ Gmᵢmⱼ×(x̲ⱼ-x̲ᵢ)×(1/|(x̲ⱼ-x̲ᵢ)|³)
	ẋ̲ᵢ = dx̲ᵢ/dt = p̲/mᵢ
which will be numerically integrated with the Runge-Kutta 4 Method to approximate the bodies trajectories.

Author: Melanie Heinrich, based on the code '/Ingolstadt.jl/tree/main/src/Development/NBodies/NBodies.jl' 
published on github by Niall Palfreyman, downloaded: 26/5/2022. 
"""
module NBodies

using GLMakie
using LinearAlgebra
#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	NBody
An NBody system capable of containing multiple (N) bodies that gravitationally interact
with each other.
"""
mutable struct NBody
	N							# Number of bodies
	nsteps						# Duration of simulation
	dt							# Timestep resolution
	G							# Gravitational constant
	x0							# Initial positions of bodies
	p0							# Initial momenta of bodies
	m							# Masses of bodies
	t							# Current time of system
	x							# Current position of system
	p							# Current momentum of system
	Ekin						# Current kinetic energy of system
	Epot						# Current potential energy of system
	Etot						# Current total energy of system		
	cp							# Current value of chaos test translation variable cp 
	cq							# Current value of chaos test translation variable cq

	"Construct a new NBody"
	function NBody( T=40, resolution=20000, G=1)
		# Initialise all fields of the decoding apparatus:
		new(
			0,					# Initially no bodies in the system
			resolution,			# Duration
			T/resolution,		# Timestep resolution
			G,					# Gravitational constant
			[],					# Initial positions
			[],					# Initial momenta
			[],					# Masses of bodies
			0.0,				# Start time at zero
			[],					# Current position empty
			[],					# Current momentum empty
			[],					# Current Ekin
			[],					# Current Epot
			[],					# Current Etot
			[],					# Current cp
			[]					# current cq
		)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	addbody!( nbody::NBody, x0::Vector{Float64}, p0::Vector{Float64}, m::Float64=1)
Add to the nbody system a new body with initial position and momentum x0, p0, and with mass m.
"""
function addbody!( nbody::NBody, x0::Vector{Float64}, p0::Vector{Float64}, m::Float64=1.0)
	push!( nbody.x0, x0)
	push!( nbody.x, deepcopy(x0))					# deepcopy() to copy by value and not by reference 
	push!( nbody.p0, p0)						
	push!( nbody.p, deepcopy(p0))
	push!( nbody.m, m)
	nbody.N  += 1
end

#-----------------------------------------------------------------------------------------
"""
	relpos( locations::Vector)
Internal utility function: Calculate matrix of antisymmetric relative positions of the Bodies.
"""
function relpos( locations::Vector)
	locPerBody = repeat(locations,1,length(locations))		# Set up locations alongside each other
	locPerBody - permutedims(locPerBody)					# Antisymmetric relative positions
end

#-----------------------------------------------------------------------------------------
"""
	massproduct( masses::Vector, G=1.0)
Internal utility function: Calculates mass product with given mass and gravitational constant.
"""
function massproduct( masses::Vector, G=1.0)
	G * masses * masses' 								# Gravitational constant * mᵢⱼ Matrix
end

#-----------------------------------------------------------------------------------------
"""
	forceOnSources( locations::Vector, gmm::Matrix{Float64})
Internal utility function: Calculate the force on each source with given massproduct at
given locations.
"""	 
function forceOnSources( locations::Vector, gmm::Matrix{Float64}) 
	relloc = relpos(locations)								
	invCube = abs.(relloc'.*relloc) .^ (3/2)				# (1/|(x̲ⱼ-x̲ᵢ)|³)
	for i in 1:length(locations) invCube[i,i] = 1 end		# Prevent zero-division along diagonal
	forceFROMSources = -gmm .* relloc ./ invCube			# Force contribution FROM each source
	sum( forceFROMSources, dims=2)							# Sum all forces ON each source
end

#-----------------------------------------------------------------------------------------
"""
	rkstep!( nbody::NBody)
To approximate solution curve of ṗ̲ᵢ & ẋ̲ᵢ of given n-body system, the Runge-Kutta 4th order method 
estimates the next step x̲ₙ₊₁ & p̲ₙ₊₁ using current x̲ₙ & p̲ₙ and weighted sum of 4 slopes.
Over a duration divided into res timesteps.
"""
function rkstep!( nb::NBody)
	dt = nb.dt							
	dt2 = nb.dt/2										# Half time step
	gmm = massproduct(nb.m)								# Gmᵢmⱼ for force calculation

	k1x = nb.p ./ nb.m									# k1: Slope at beginning of timestep 							
	k1p = forceOnSources(nb.x, gmm)[:]

	k2x = (nb.p .+ (dt2*k1p)) ./ nb.m					# k2: Estimated slope at midpoint using k1
	k2p = forceOnSources((nb.x + (dt2*k1x)), gmm)[:]	 

	k3x = (nb.p .+ (dt2*k2p)) ./ nb.m					# k3: Estimated slope at midpoint using k2
	k3p = forceOnSources((nb.x + (dt2*k2x)), gmm)[:]	

	k4x = (nb.p .+ (dt*k3p)) ./ nb.m					# k4: Estimated slope at endpoint using k3
	k4p = forceOnSources((nb.x + (dt*k3x)), gmm)[:]	

	# Weighted sum of slopes for final estimates:
	nb.x = nb.x + dt * (1/6) * (k1x +2k2x +2k3x +k4x)	# x̲ₙ₊₁
	nb.p = nb.p + dt * (1/6) * (k1p +2k2p +2k3p +k4p)	# p̲ₙ₊₁

end

#-----------------------------------------------------------------------------------------
"""
	getenergy!( nbody::NBody)
Calculate current kinetic-, potential- and total-energy with the given positions and 
momenta of the bodies in the system. With following equations:
	Eₖᵢₙ = ∑ᴺᵢ₌₁ (1/2)×mᵢ×|v̲ᵢ|² = ∑ᴺᵢ₌₁ (1/2)×|p̲ᵢ|²×(1/mᵢ)
	Eₚₒₜ = -∑ᴺᵢ₌₁ ∑ᴺⱼ₌ᵢ₊₁ Gmᵢmⱼ×|(x̲ⱼ-x̲ᵢ)|
	Eₜₒₜ = Eₖᵢₙ + Eₚₒₜ
"""
function getenergy!( nb::NBody)	
	# Kinetic Energy
	pmag = Vector{Float64}(undef,length(nb.p))			# Holds magnitude of current momentum for each body
	for i in 1:length(nb.p) 							
		pmag[i] = sqrt.(sum(nb.p[i].^2))				# |p̲ᵢ| = √(p(x)ᵢ² + p(y)ᵢ²)
	end
	nb.Ekin = sum((1/2)*((pmag.^2) ./ nb.m))			# Eₖᵢₙ = ∑ᴺᵢ₌₁ (1/2)×|p̲ᵢ|²×(1/mᵢ) 

	# Potential Energy
	gmm = massproduct(nb.m)								# Gmᵢmⱼ 
	relloc = relpos(nb.x)
	xmag_ji = sqrt.(abs.(relloc'.*relloc))				# |(x̲ⱼ-x̲ᵢ)|
	for i in 1:length(nb.x) xmag_ji[i,i] = 1 end		# Prevent zero-division along diagonal
	epot_matrix = gmm ./ xmag_ji						# Eₚₒₜ = -∑ᴺᵢ₌₁ ∑ᴺⱼ₌ᵢ₊₁ Gmᵢmⱼ×|(x̲ⱼ-x̲ᵢ)|
	nb.Epot = -sum(triu(epot_matrix))					# Sum over upper triangle to count each interaction only once (j=i+1)
	
	# Total Energy
	nb.Etot = nb.Ekin + nb.Epot

end

#-----------------------------------------------------------------------------------------
"""
	pq_chaostest!( nbody::NBody, n::Int)
Calculate the translation variables cp and cq of the 0-1 test for chaos with the current positions.
A one dimensional time series 𝜑(n)(here current x̲ₙ) undergoes Euclidean extension deriving a pair of
extended coordinates cp(n) & cq(n), parameterized by an angle c. With c being a fixed number c ∈ [0,2π].
Plotting cp vs cq derived from the current x-coordinate of the bodies, will give indication of the systems dynamics:
	Regular underlying dynamics of system: cp vs cq will show bounded dynamics.
	Chaotic underlying dynamics of system: cp vs cq will show diffusive dynamics.
"""
function pq_chaostest!(nb::NBody, n::Int)
	c = 1		# c is 1 as cp&cq are calculated with cosine and sine in degrees
	if n == 1
		nb.cp = nb.x .* cosd(c)				# cp(1) = 𝜑(n)×cos(c)
		nb.cq = nb.x .* sind(c)				# cq(1) = 𝜑(n)×sin(c)
	else
		nb.cp = nb.cp + nb.x .* cosd(c*n)	# cp(n+1) = p(n) + 𝜑(n)×cos(c×n)
		nb.cq = nb.cq + nb.x .* sind(c*n)	# cq(n+1) = q(n) + 𝜑(n)×sin(c×n)
	end

end

#-----------------------------------------------------------------------------------------
"""
	simulate( nb::NBody)
Run a simulation of the given NBody system over duration nb.nsteps * nb.dt.
"""
function simulate( nb::NBody)
	t = 0:nb.dt:nb.nsteps*nb.dt
	x = Vector{typeof(nb.x0)}(undef,nb.nsteps+1)
	p = Vector{typeof(nb.p0)}(undef,nb.nsteps+1) 
	Ekin = Vector{Float64}(undef,nb.nsteps+1)	
	Epot = Vector{Float64}(undef,nb.nsteps+1)
	Etot = Vector{Float64}(undef,nb.nsteps+1)
	cp = Vector{typeof(nb.x0)}(undef,nb.nsteps+1)
	cq = Vector{typeof(nb.x0)}(undef,nb.nsteps+1)

	# Initialisation:
	x[1] = nb.x0
	p[1] = nb.p0

	# Simulation:
	for n = 1:nb.nsteps
		# Calculate energy of system with current position and momentum
		getenergy!(nb)
		Ekin[n] = nb.Ekin
		Epot[n] = nb.Epot
		Etot[n] = nb.Etot
	
		# Calculate translation variables for current postion
		if n == 1					
			pq_chaostest!(nb, n)		# Initialisation of cp & cq
			cp[1] = nb.cp			
			cq[1] = nb.cq
		end
		pq_chaostest!(nb, n+1)
		cp[n+1] = nb.cp
		cq[n+1] = nb.cq

		# Calculate next step of position and momentum using Runge-Kutta 4
		rkstep!(nb)
		x[n+1] = nb.x
		p[n+1] = nb.p
	end

	(t, x, p, Ekin, Epot, Etot, cp, cq)
end

#-----------------------------------------------------------------------------------------
"""
	sim_anim_plot( nbody::NBody)
Simulates given 4-body system. Animates and plots: 
	- Current locations to animate movement of bodies along their trajectories
	- Current energy of the system(Eₖᵢₙ, Eₚₒₜ, Eₜₒₜ) along the plotted graph over all time steps
"""
function sim_anim_plot( nb::NBody)
	# Run the simulation:
	(t, x, p, Ekin, Epot, Etot, cp, cq) = simulate(nb)

	# Observables for animation
	x_curr = Observable(map(body-> body[1],x[1]))
	y_curr = Observable(map(body-> body[2],x[1]))
	t_curr = Observable(string("t = ", round(t[1], digits=2)))						
	Ekin_curr = Observable([Ekin[1]])
	Epot_curr = Observable([Epot[1]])
	steps_plot = [0:1:nb.nsteps]												
	step_plot_curr = Observable([map(point->point, steps_plot[1][1])])		# To plot against energy observables

	# Plotting and animation of the body system:
	fig = Figure(resolution=(1200, 700))
	ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "N-body 2D Motion on their Trajectories")
	limits!(ax1, -1.2, 2.2, -1.2, 2.2)
	# Trajectories plotting
	lw1 = 1.5						# Linewidth
	lines!(ax1, map(el->el[1][1],x), map(el->el[1][2],x), linewidth = lw1, color=:plum)					# Body 1
	lines!(ax1, map(el->el[2][1],x), map(el->el[2][2],x), linewidth = lw1, color=:lightskyblue3)		# Body 2
	lines!(ax1, map(el->el[3][1],x), map(el->el[3][2],x), linewidth = lw1, color=:peachpuff3)			# Body 3
	lines!(ax1, map(el->el[4][1],x), map(el->el[4][2],x), linewidth = lw1, color=:darkseagreen2)		# Body 4
	# Animation with current values
	scatter!(ax1, x_curr, y_curr, markersize=(20), color=[:plum, :lightskyblue3, :peachpuff3, :darkseagreen2])
	text!( t_curr, position = (-0.5, 1.5),fontsize=18, align=(:center, :top))

	# Energy diagram graph plotting
	ax2 = Axis(fig[1,2], xlabel = "Duration in n timesteps", ylabel = "Energy in Joule", title = "Energy of the System")
	limits!(ax2, 0, nb.nsteps, -60, 60)	
	lw2 = 2.5						# Linewidth
	lines!(ax2, Ekin, label="Kinetic energy", linewidth = lw2, color=:darksalmon)			# Ekin
	lines!(ax2, Epot, label="Potential energy", linewidth = lw2, color=:darkseagreen)		# Epot
	lines!(ax2, Etot, label="Total energy", linewidth = lw2, color=:grey33)					# Etot
	axislegend(ax2, labelsize=13, position=:lt)
	# Energy diagram animation with current values
	scatter!(ax2, step_plot_curr, Ekin_curr,  markersize=(13), color=[:darksalmon])
	scatter!(ax2, step_plot_curr, Epot_curr, markersize=(13), color=[:darkseagreen])

	display(fig)
	#Run the animation:
	for i in 1:length(t)-1
		x_curr[] = map(body->body[1],x[i])
		y_curr[] = map(body->body[2],x[i])
		t_curr[] = string("t = ", round(t[i]; digits=1))
		Ekin_curr[] = [Ekin[i]]	
		Epot_curr[] = [Epot[i]]	
		step_plot_curr[] = [i]
		sleep(1e-12)
	end

end
#-----------------------------------------------------------------------------------------
"""
	sim_plot( nbody::NBody, fig::Makie.Figure, i::Int)
Simulates given (2 to 4)-body system and places plots on given Figure at position i: 
	- Trajectories of bodies over all time steps 
	- Translation variables cp versus cq from the 0-1 test for chaos
"""
function sim_plot( nb::NBody, fig::Makie.Figure, i::Int)
	# Run the simulation:
	(t, x, p, Ekin, Epot, Etot, cp, cq) = simulate(nb)

	# Trajectories plotting of the body system:
	lw1 = 1.5						# Linewidth
	ax1 = Axis(fig[1,i], xlabel = "x", ylabel = "y", title = "N-body 2D Trajectories")							# Upper plot at position i
	limits!(ax1, -6, 6, -4, 8)
	lines!(ax1, map(el->el[1][1],x),map(el->el[1][2],x),linewidth = lw1, color=:plum)								# Body 1
	lines!(ax1, map(el->el[2][1],x), map(el->el[2][2],x),linewidth = lw1, color=:lightskyblue3)						# Body 2
	if nb.N >= 3 lines!(ax1, map(el->el[3][1],x), map(el->el[3][2],x),linewidth = lw1, color=:peachpuff3)	end		# Body 3
	if nb.N == 4 lines!(ax1, map(el->el[4][1],x), map(el->el[4][2],x),linewidth = lw1, color=:darkseagreen3) end	# Body 4

	# 0-1 test for chaos: Plotting cp versus cq only derived of the x-coordinates for each body 
	ax2 = Axis(fig[2, i],xlabel = "cp", ylabel = "cq", title = "cp vs. cq from the 0-1 Test for Chaos")			# Lower plot at position i
	limits!(ax2, -1000, 1000, -1000, 1000)
	scatter!(ax2, map(el->el[1][1],cp), map(el->el[1][1],cq), markersize=(1), color=:plum)							
	scatter!(ax2, map(el->el[2][1],cp), map(el->el[2][1],cq), markersize=(1), color=:lightskyblue3)
	if nb.N >= 3 scatter!(ax2, map(el->el[3][1],cp), map(el->el[3][1],cq), markersize=(1), color=:peachpuff3) end		
	if nb.N == 4 scatter!(ax2, map(el->el[4][1],cp), map(el->el[4][1],cq), markersize=(1), color=:darkseagreen3) end

end

end		# of NBodies