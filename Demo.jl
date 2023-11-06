"""
Main which uses the module NBodies to demonstrate the simulation of a 4 body-system and show the influence 
of only slightly differing masses of the bodies on the periodic or chaotic nature of the system.
"""
module Demo

include("NBodies.jl")
using .NBodies
using GLMakie

"""
	demo()
Demonstrate simulation of 3 different n-body systems.
"""
function demo()
	# Initialise a 4-body system which is periodic for awhile:
	nb1 = NBodies.NBody( 20,4000, 1)
	NBodies.addbody!( nb1, [0.00,0.00],[1.00,0.00],1.00)
	NBodies.addbody!( nb1, [0.00,1.00],[0.00,-1.00],1.00)
	NBodies.addbody!( nb1, [1.00,1.00], [-1.00,0.00], 1.00)
	NBodies.addbody!( nb1, [1.00,0.00], [0.00,1.00], 1.00)

	displaytext1 = """This demo function demonstrates the simulation of a 4 body-system and shows 
	the influence of only slightly differing the masses of the bodies on the periodic or chaotic nature of the system.

	First will be the simulation of a simple 4-body system with unit mass which is periodic for a very short period.
	At some point two bodies are about to collide which changes the total energy displayed for the system.

	"""
	displaytext2 = """
	The system should not gain or loose energy, the problem is the numerical integration of the equations of motion.  
	As the two point masses are about to collide, their relative positions come extremly close to zero. 
	With F̲ᵢⱼ  = ∑ Gmᵢmⱼ×(x̲ⱼ-x̲ᵢ)×(1/|(x̲ⱼ-x̲ᵢ)|³), in the next step calculated with Runge-Kutta Method the acceleration skyrockets.
	Resulting in calculating extremly high magnitudes for current potential and kinetic energy and a change in total energy.
	From that moment on, the energy calculated for the system is wrong.

	"""
	println(displaytext1, ">>Press Enter in Julia Repl to load animation & energy diagram.<<")
	# Wait for user input
	readline()
	println(displaytext2)
	# Run simulation, animate & plot trajectories and energy of the 4-Body-System nb1
	NBodies.sim_anim_plot(nb1)
	displaytext3 = """The system is very unstable. There are many 4-body systems such as stellar systems which are periodic, 
	the big difference is their bodies mass differs from each other.
	
	3 different n-body systems will be simulated for the same time & duration to show the influence of differing masses of the bodies.
	Beneath their plotted trajectories are plots to test for the chaotic or periodic nature of the system,
	cp plottet versus cq are bounded if the underlying dynamics is regular, i.e. periodic or quasiperiodic.
		Left: Same system as previously simulated.(all unit mass)
		Middle: Same parameters as the system on the left, but with slightly differing masses.
			Masses of bodies: purple:0.995, green:0.91, blue:1.0002, beige:1.000151  
		Right: 3-body problem solution which is periodic.

	"""
	println(displaytext3, ">>Press Enter in Julia Repl to show plots.<<")
	readline()
	# Initialising nb1 again for same duration as nb2 & nb3
	nb1 = NBodies.NBody( 80, 8000, 1)
	NBodies.addbody!( nb1, [0.00,0.00], [1.00,0.00], 1.00)
	NBodies.addbody!( nb1, [0.00,1.00], [0.00,-1.00], 1.00)
	NBodies.addbody!( nb1, [1.00,1.00], [-1.00,0.00], 1.00)
	NBodies.addbody!( nb1, [1.00,0.00], [0.00,1.00], 1.00)
	# Initialise a 4-body system which is more of a periodic nature than nb1:
	nb2 = NBodies.NBody( 80, 8000, 1)
	NBodies.addbody!( nb2, [0.00,0.00], [1.00,0.00], 0.995)		    # Purple
	NBodies.addbody!( nb2, [0.00,1.00], [0.00,-1.00], 1.0002)	    # Blue
	NBodies.addbody!( nb2, [1.00,1.00], [-1.00,0.00], 1.000151)	    # Beige
	NBodies.addbody!( nb2, [1.00,0.00], [0.00,1.00], 0.91)	        # Green
	# Initialise a periodic and stable 3-body system
	nb3 = NBodies.NBody( 80, 8000, 1)
	NBodies.addbody!( nb3, [1.666163752077218-1,-1.081921852656887+1], 	[0.841202975403070,0.029746212757039],	1.00)		
	NBodies.addbody!( nb3, [0.974807336315507-1,-0.545551424117481+1],	[0.142642469612081,-0.492315648524683], 1.00)		
	NBodies.addbody!( nb3, [0.896986706257760-1,-1.765806200083609+1],	[-0.983845445011510,0.462569435774018], 1.00)	
	
	figure = Figure(resolution=(1400, 900))
	# Run simulation, plot trajectories and chaos plots of 3 different n-body systems
	NBodies.sim_plot(nb1, figure, 1)
	NBodies.sim_plot(nb2, figure, 2)
	NBodies.sim_plot(nb3, figure, 3)
	display(figure)
end

end # ofDemo
