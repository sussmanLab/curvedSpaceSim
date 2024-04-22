#Test Explainer
Each of the tests is essentially run through command line inputs, 
parsed via TCLAP. They're all variations on the core sample code 
in curvedSpaceSimulation.cpp. For each of the time vs. particle
number and time vs. mesh complexity tests, we strictly 
utilized harmonically repulsive forces for flexible numbers of 
iterations at each value tested. 

##Cost vs. N 
Once the code compiles and you are able to cmake -> make 
successfully, one should be able to run the cost vs. 
N test totally from the command line, with minor changes
to the code per the user's objective. 

The command line arguments that can be specified for
this test are the number of iterations for each 
particle number, the number of particle numbers
to try (referred to as doublings, as we used
multiples of two for our tests), the mesh one 
wants to use, the interaction range for the particle
interactions (especially relevant for submeshed tests,
as this sets the submesh length scale), the timestep
size, and whether to use a reproducible RNG.  

If one were to specify all of those via TCLAP with the settings
we used to make figure 3, the command line input
would be: 
```./t_v_n.out -i 5 -n 15 -m "../exampleMeshes/torus_isotropic_remesh.off" -a 1 -t .01 -r 1```

This input will run both all-to-all and submeshed tests
using the interaction range specified by -a. The all-to-all
tests take quite a long time (and in fact, they are the limiting
factor on the particle number scale we can use for these tests), 
so if one only wishes to reproduce the scaling behavior
for submeshed tests, please comment out the all-to-all test
at the designated comment markers. One file will be produced
with the particle numbers and paired average simulation timestep 
costs. 

##Cost vs. Complexity


##Crystallization Simulations

