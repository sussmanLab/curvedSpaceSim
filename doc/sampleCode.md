# Code snippets

These minimal samples show how to run default-initialized simulations. Please see the provided .cpp
files in the main directory

# curvedSpaceSimulation.cpp

Basic simulations on (by default) a torus.

# multirankSimulation.cpp

Basic MPI simulations on (by default) a torus. Called via things like
$mpirun -n XX ./multirankSimulation.out [command-line arguments here]
where XX is the number of processes you'd like to run on.

## meshTesting.cpp and flatSpaceSimulation.cpp

Two files that exist largely for debuging and testing code, either by checking 
specific mesh-based functionality or by running simulation algorithms in flat (R^3)
spatial domains rather than on curved surfaces.
