# spatial_decision_making
Here is the replication of the neural model of "The geometry of decision-making in individuals and collectives" Vivek H. Sridhar, Liang Li, Dan Gorbonos, Máté Nagy, Bianca R. Schell, Timothy Sorochkin, Nir S. Gov and Iain D. Couzin https://doi.org/10.1073/pnas.1811964115

## Background

Collective decision-making occurs at various biological scales: both groups of individuals and collections of neurons integrate information to make decisions. In addition, the geometry of a decision-making problem could be key to understand the mechanisms underlying movement directions. 
This very cool paper uses both theory and experiments to investigate how geometry could influence decision making at mulitple biological scales, both collections of individuals and collections of neurons. Specifically, when the angular difference between attracting options is below a critical threshold, collectives average the possible directions, but once the angular difference exceeds the critical threshold, the collective breaks symmetry and makes a consensus decision towards one of the options. The incresed sensitivity near cricitcality might also provide decision-making advantages, which despite often suggested has never been definetly proved. This cool paper goes very well in that direction! I have replicated the neural decision-making model, which implements a Hopfiel network, a simple type of fully connnected binary network which have been extensively used in pyshics, molecular biology, and neurosciences. 

## Code
The neural network has a spatial position and has to move towards two possible targets. The network is composed of "spins" (not called "neurons" to avoid confusing the Hopfield network with a mechanistic rappresentation of an animal brain) which encode the preferred direction to one of two targets and can be either turned on or off. The weights of the edges between the spins that are turned on sum up and give the energy function of the pattern of spin activation (called Hamiltonian). The energy decreases if the activated neurons encode similar prefered directions. Every time step, the energy function is minimised thorugh an optimization technique called Metropolis-Hasting algorithm by turning on or off the spins. The network than "moves" a tiny bit according to the average direction encoded by its activated spins, and the process is repeated untill the network reaches one of the targets. Since it was written in R, the functions handling spatial calculations (like finding heading and an angular difference between vectors) have been vectorized for computational efficency. Also calculating the energy function relies on R vectorization to avoid for loops. Nevertheless, the simulation is quite slow, and only a few replicates of the movement trajectory of the network are shown.
