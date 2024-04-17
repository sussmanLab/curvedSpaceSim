# How to Write a Simulation

This file aims to explain how one would write a custom simulation file 
by highlighting the key parts of the example code. To do so, it will reference
curvedSpaceSimulation.cpp and explain why the file is structured in the way 
it is. 

##Structure of a Simulation Executable

We can break a simulation executable file down into its core parts:
1. Command line arguments. These can be anything from a control on whether submeshing is done
   to temperature controls. To use a command line variable, we first create a ValueArg 
   of the desired type, then assign the value to a variable via .getValue(). After this, 
   we can use the variable normally as if it had been defined within our main function. 
2. Space. This can be e.g. flat Cartesian space or a space defined by a mesh. If using a mesh
   space, one must specify an input file and load it using loadMeshFromFile. Under the hood, 
   this is just a call to CGAL::PMP's read polygon mesh, which takes .off, .obj, .stl, .ply, 
   .ts, and .vtp files. Additionally, the space contains the means by which to measure distance,
   so it is where we decide whether to submesh or not (use a cutoff to avoid large sequence trees)
   when calculating distance. If submeshing, one must also define a cellList neighbor structure for
   the simulation to use; this relies on the extent of the underlying mesh. 
3. Configuration. This is where one defines particle positions. Configurations are tied to
   meshSpaces (as positions exist within associated spaces) via setSpace. One can establish a 
   set of positions in two ways: first, one can make a noiseSource object and then use it 
   setRandomParticlePositions, which is a great way to perform preliminary investigations of 
   a space. The particle positions are set as random barycentric weights in a random face of the
   mesh (so randomly defined meshPosition objects).  
4. Force. This is how one defines, at present, interparticle forces that are used in displacement
   calculations. The specific details of force instantiation depend on the force being used. 
   a force is also one of the easiest things to add yourself -- you need only define energy and 
   force values, as can be seen in, for example, harmonicRepulsion.cpp and harmonicRepulsion.h.
   Importantly, a simulation can be defined with multiple forces. If you are experiencing suspiciously
   slow timesteps, ensure that you are not evaluating too many forces at once; for example, if you only wish to use
   one force but have previously defined another, you should call the simulation's
   clearForceComputers() before adding in the desired force function. 
5. Updater. The updater defines how to update particle positions. It will typically call a force function,
   but might also involve some other rules about how displacements are incorporated. Familiar examples 
   from other simulation techniques are gradient descent, velocity Verlet (both ofwhich we implement here),
   or FIRE (which we might implement in the future). This is also where one would likely add unique update
   rules, like Vicsek updates.  
6. Simulation object. The simulation is like a wrapper; it stores the configuration, space, force,
   and updater. All of those things, being tied together by the simulation object, can access 
   the key pieces from each other, but also be individually modified. 
7. Main loop. This is where we actually use the simulation object that we've created. Usually, 
   one will save the trajectory into a netCDF data file to be loaded later in the environment of
   choice. It's best not to save too frequently -- trajectory files can be very large, even for 
   relatively small particle numbers. The example in the curvedSpaceSimulation.cpp uses a for loop 
   that terminates upon reaching a certain maximum number of iterations, but one might also consider,
   for example, a while loop that terminates upon meeting some condition (e.g. a very low value of 
   the forceNorm, which we get out to assess simulation status in curvedSpaceSimulation.cpp). 


##Walkthrough
