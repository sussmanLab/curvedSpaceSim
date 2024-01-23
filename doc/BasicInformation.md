# Basic overview of the project

This project is, in many ways, a very standard implementation of molecular-dynamics-
like particle-based simulations. The primary change is that all distances are computed
as geodesic distances on curved surfaces that have been approximated by a triangulated
mesh. Along with that, all particle displacements occur along discrete geodesic paths,
and things like velocity vectors are parallel transported as well. This is made possible
by interfacing the code with an implementation of an exact discrete geodesic algorithm,
and it is made computational feasible by implementing some aggresive "submeshing" of
the full surface so that finding geodesic distances can be done at a reasonable cost.

## Classes and data structures of interest

### Primary data structure

At the moment the primary data structure is a triple of doubles together with an integer.
These typically encode the index of the face a particle is on and the barycentric 
coordinates of the particle in that face, although in some data structures (e.g.,
in the cellListNeighborStructure) the three doubles are in fact the R^3 position 
of the particle and the int is meaningless.

### src/models directory

This directory contains both "spaces" and "models." Spaces are classes that know how
to compute distances (geodesic or otherwise), displace positions, and transport 
vectors. Models are classes that have particle-like data as members, can access a
space to do the operations above, and might have access to a structure that helps find
nearby neighboring particles that might be interacting with a given particle.
simpleModels are used most of the time, but mpiModels have implemented the simplest
possible version of spreading the workload across multiple processes, splitting
particle indices among the available threads by sharing global information about all
particle positions at each step.


### src/updaters directory

Contains the classes that implement various equations of motion that the simulation dynamics can
evolve according to.

### src/forces directory

Contains classes for computing interparticle forces. At the moment implementing new
pairwise interactions is quite easy; new many-body forces will require more work.

### src/simulation directory

"simulations" (basic or mpi) are classes that tie together a model (with its space),
forces, and updaters so that a standard "timestep" can be easily called with one line.

### src/utilities directory

Contains assorted utilities (access to various mesh functions, structures that
accelerate the process of finding neighbors, submeshing helpers, random-number
generators, etc).

### src/databases  directory

For convenience, Contains structures that may be useful for saving simulation
trajectories. ost of the databases are currently either simple txt files or netCDF.

### inc directory

Contains a few standard headers (defining point data types, useful debugging commands,
etc).

### exampleMeshes directory

For convenience, we include a small handfull of triangulated mesh surfaces to play with.
