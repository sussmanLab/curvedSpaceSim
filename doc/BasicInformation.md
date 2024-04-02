# Basic overview of the project

This project is, in many ways, a very standard implementation of molecular-dynamics-
like particle-based simulations. The primary change is that all distances are computed
as geodesic distances on curved surfaces that have been approximated by a triangulated
mesh. Along with that, all particle displacements occur along discrete geodesic paths,
and things like velocity vectors are parallel transported as well. This is made possible
by interfacing the code with an implementation of an exact discrete geodesic algorithm,
and it is made computationally feasible by implementing an aggresive "submeshing" of
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

### src/simulation directory

"simulations" (basic or mpi) are classes that tie together a model (with its space),
forces, and updaters so that a standard "timestep" can be easily called with one line.
The simulation object is like a wrapper to store important details about the 
simulation you want to run such that they're easily accessible from within executable 
files. 

### src/updaters directory

Contains the classes that implement various equations of motion that the simulation dynamics can
evolve according to. Updaters are child classes of a basic updater, named baseUpdater. 
The only requisite component of an updater is a performUpdate() function, which indicates
how that updater is supposed to update particle positions. Typically, this involves the
computation of "forces," which are dependent on the force currently loaded into the 
simulation object you're using. Following force computation, the updater will call
a displacement function -- the exact form of which will depend on the model space 
stored in the simulation object. For example, if using a triangulatedMeshSpace, 
the displacement function will involve bending the path around the mesh. 

Updaters can also include more complex update rules. For example, the Verlet updater
incorporates velocity contributions to updates as well as force evaluations. In general,
if you want to change the physical model that you're simulating, you will want to change
either an updater or a force. For larger structural changes to the model, you will likely
aim to change an updater. 

### src/forces directory

Contains classes for computing interparticle forces. At the moment implementing new
pairwise interactions is quite easy; new many-body forces will require more work.

Pairwise interactions are taken as implicit assumptions in, for example, the gaussianRepulsion
and harmonicRepulsion force calculators, which take only distances and separation vectors
as inputs. A separation vector indicates a vector along the tangent to the path between
two points. Current sign conventions always point from source to target on both ends
of the straight-line or geodesic path -- forces on source particles therefore always
include a negative sign.    

### src/utilities directory

Contains assorted utilities (access to various mesh functions, structures that
accelerate the process of finding neighbors, submeshing helpers, random-number
generators, etc).

For a description of the individual functionalities: 
- Neighbor structures: 
	- baseNeighborStructure: The basic way of organizing particle neighbors. 
	  It defaults to all-to-all lists of neighbors, but includes important inheritable
	  features for non-all-to-all approaches. For example, a way to set the interaction 
	  range, initialization of the structure, and a means by which to construct the list
	  of candidate neighbors (which is constructed via Euclidean distance and trimmed by 
	  geodesic distance).  
	- cellListNeighborStructure: Technically also an inheritable class, although this
	  one doesn't assume all-to-all interactions -- instead, it assumes we will organize
	  space via a cell list (a list of space partitions were neighbors might exist). 
	  To accomplish this, cellListNeighborStructure.cpp marks out the minimum and maximum
	  positions in space when it is instantiated, establishes a grid size using those 
	  spatial boundaries and a grid size set by the interaction range, and finds the 
	  list of candidate neighbors. 
- Cell list structure, hyperRectangularCellList: partitions space into a grid of 
  rectangles with side lengths defined by the extent of space covered by an underlying
  triangulatedMeshSpace (or free space, if you so desire). The indexer is used here to 
  assign indices to the cells. After space has been partitioned, we can use its partioning
  to manage neighbors; getCellNeighbors determines which cells are the index of the 
  cellIndex we care about, while sort arranges meshPositions into their cells.  
- functionUtilities: A variety of simple utilities used by functions:
	- vectorMagnitude and normalize, to determine and normalize vector sizes
	- rotateAboutAxis, which rotates a point about a given axis defined by two other points
	- angleBetween, which determines the angle between two vectors
	- wrap, fitting integers into domains
- meshUtilities: The bulk of mesh manipulation and analysis tools used by other functions.
  Most have self-explanatory names.
	- meshPositionToFaceLocation turns a meshPosition in our code into a Face location
	  type for certain CGAL functions, while getVertex (index/position) From Face yields
	  either the original index of the face's vertices in the mesh or their XYZ positions 
	- meanEdgeLength and meanTriangleArea give the named quantities from the mesh's polygons
	- Here is where the intersection checker is (necessary for displacements) -- 
	  it looks complicated because of the moving
	  pieces, but it's genuinely just a barycentric line intersection checker which is 
	  complex because of the several different edge-related quantities that must be 
	  calculated. See https://math.stackexchange.com/questions/3995636/point-of-intersection-of-two-lines-in-barycentric-coordinate-system. 
	- Finally, meshUtilites is where the functionality to compute distances and tangents
	  (computePathDistancesAndTangents) is located. It's later used in force functions and 
	  so forth; this is a call to the CGAL Xin and Wang algorithm with a given mesh 
	  that can either be the original mesh of a submesh of it. 
- submesher: Uses a source point and a set of targets (the confirmed neighbors of the source point)
	     to create a submesh. To do this, we do a depth-first search of all the faces 
	     in the interacting region that might contain a target (added to the 
	     explorationStack). In doing so, we're sure to keep the region connected so that 
	     there are no holes which might distend geodesics. 
	     After we've finished visiting faces that might contain target points, we 
	     pass their indices up to a function which takes the information of those faces
	     in the original mesh and generates a new mesh, specifically using vertex and 
	     edge information.  
	     
- noiseSource: An easy way of getting random numbers. The most useful function for simple
  simulations is getRandomBarycentricSet, which gives a trio of doubles that act as random
  barycentric coordinates; pair those with a random int, and you have a random set of 
  barycentric coordinates on a random face, as in, the means by which to get uniform 
  samples across a mesh. 
- indexer: Helper header file for assigning indices, especially to cells in 
  hyperRectangularCellList. Both 2-D and 3-D indexing is supported. Most importantly,
  this lets us convert a 2- or 3-D array of indices to a much faster 
  one-dimensional equivalent.
- profiler: Specialty timer for profiling code functionality. One can define a profiler
  as a class profiler with this header and then wrap a function with its start and end 
  functions; these call on c++'s chrono functionality to get current times and store 
  differences/number of function calls/etc. 
   


### src/databases  directory

For convenience, contains structures that may be useful for saving simulation
trajectories. ost of the databases are currently either simple txt files or netCDF.

### inc directory

Contains a few standard headers (defining point data types, useful debugging commands,
etc).

### exampleMeshes directory

For convenience, we include a small handfull of triangulated mesh surfaces to play with.
