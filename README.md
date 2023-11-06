# N-body-simulation
### In this module a system of N interacting bodies moving in 3 dimensions over time will be simulated. The trajectories of the bodies are approximated by numerically integrating the equations of motion with the Runge-Kutta Method.

## N-body equations of motion are defined by: 
###   F̲ᵢⱼ = mᵢ×(d²x̲ᵢ/dt²) = ∑ Gmᵢmⱼ×(x̲ⱼ-x̲ᵢ)×(1/|(x̲ⱼ-x̲ᵢ)|³), where the right side is summed up over N with j=1 & j≠i. 
These equations describe the motions of N mass points mᵢ, 
moving under the influence of their mutual attracting force given by Newton's law of gravitation.
with p̲ᵢ = mᵢ×(dx̲ᵢ/dt) follows a set of first order differential equations:
###   ṗ̲ᵢ = dp̲ᵢ/dt = ∑ Gmᵢmⱼ×(x̲ⱼ-x̲ᵢ)×(1/|(x̲ⱼ-x̲ᵢ)|³)
###   ẋ̲ᵢ = dx̲ᵢ/dt = p̲/mᵢ
which will be numerically integrated with the Runge-Kutta 4 Method to approximate the bodies trajectories.

Author: Melanie Heinrich, based on the code '/Ingolstadt.jl/tree/main/src/Development/NBodies/NBodies.jl' 
published on github by Niall Palfreyman, downloaded: 26/5/2022. 
