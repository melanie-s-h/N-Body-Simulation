# N-Body-Simulation
<br>

***In this module a system of N interacting bodies moving in 3 dimensions over time will be simulated. The trajectories of the bodies are approximated by numerically integrating the equations of motion with the Runge-Kutta Method.***
<br>
<ins>
**N-body equations of motion are defined by**
</ins>
<br>
###           *F̲ᵢⱼ = mᵢ×(d²x̲ᵢ/dt²) = ∑ Gmᵢmⱼ×(x̲ⱼ-x̲ᵢ)×(1/|(x̲ⱼ-x̲ᵢ)|³)*,<br>			
*where the right side is summed up over N with j=1 & j≠i*.	
<br>These equations describe the motions of N mass points mᵢ, 
moving under the influence of their mutual attracting force given by Newton's law of gravitation.	
with **p̲ᵢ = mᵢ×(dx̲ᵢ/dt)** follows a set of first order differential equations:				
<br>

###           *ṗ̲ᵢ = dp̲ᵢ/dt = ∑ Gmᵢmⱼ×(x̲ⱼ-x̲ᵢ)×(1/|(x̲ⱼ-x̲ᵢ)|³)*
<br>

###          *ẋ̲ᵢ = dx̲ᵢ/dt = p̲/m*

which will be numerically integrated with the Runge-Kutta 4 Method to approximate the bodies trajectories.
<br>

*Author: Melanie H., based on the code '/Ingolstadt.jl/tree/main/src/Development/NBodies/NBodies.jl' 
published on github by Niall Palfreyman, downloaded: 26/5/2022.* 


<hr>

<br>

Type Demo.demo() in Julia REPL to run the simulation.<br>
The "demo()" function of "Demo.jl" file demonstrates the simulation of 3 different n-body systems while giving explanation in the REPL.
"Demo.jl" uses the main module "NBodies.jl", where the math and physics is happening.
<br>


![image](https://github.com/melanie-s-h/N-body-simulation/assets/134691659/95a147fb-c79f-4f44-8561-297ffa232b24)





https://github.com/melanie-s-h/N-body-simulation/assets/134691659/8539820c-c420-4f96-99d2-2bb675dbd990


![image](https://github.com/melanie-s-h/N-body-simulation/assets/134691659/00935c8a-a8cc-4166-8871-24251f7c5ccb)

![image](https://github.com/melanie-s-h/N-body-simulation/assets/134691659/0679b63b-a9ba-4aae-96ce-cc5b42670b91)



<hr>

