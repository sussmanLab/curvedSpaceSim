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

This walkthrough essentially explains the structure of curvedSpaceSimulations.cpp,
and will generally assume the defaults for command line arguments. 

0. The target for our simulation will be gradient descent on a torus. 
1. First, we need to define our command line arguments. Because we're doing a simple simulation 
   of particles on a torus, we need to be able to import a toroidal mesh (as a string file name).
   We also need to define particle number, set ineteraction range, allow for submeshing (or not), establish how frequently we'll
   save, and figure out how long the simulation will go. That means we need two ```ValueArg<int>``` for particle
   number and save frequency, a ```ValueArg<string>``` for the mesh name, a ```ValueArg<double>``` for interaction range, 
   and a ```ValueArg<bool>``` (or int) for turning submeshing on and off. When using bools with TCLAP, you pass to the
   program as 0 for false and 1 for true. There are other ValueArgs in the curvedSpaceSimulation.cpp file, but if
   we wanted to use those 4 in our simulation, we might initialize it as: 
> ```./curvedSpaceSimulation.out -z 1 -n 20 -a 1 -s 100 -i 1000```
2. We need to fetch all of our command line arguments, so we define variables of types that match the ValueArgs
   to store them. E.g. **int** saveFrequency = saveFrequencyArg.getValue(). 
3. Now we need to begin setting up the simulation object. A simulation requires: 
* A space
* A configuration
* One or more forces
* One or more updaters
4. First, we define the meshSpace, as the configuration (of particle positions) is dependent on it. To do so, 
   we must create a shared pointer that points the meshSpace we want to use. The shared pointer will allow other
   parts of the code to easily reference the meshSpace and what it contains. Then, we can set the associated
   mesh using ```meshSpace->loadMeshFromFile(filename)```. The arrow operator ```->``` is necessary because it 
   accesses what the shared pointer points to; so what comes before the arrow is the shared pointer, and what
   comes after is the value or function of what that pointer points to that we want to access. 
5. If we want to use restricted interaction ranges, we need to set a cellList. To instantiate a cellListNeighborStructure, 
   we have to tell the program how big the space is and how coarsely to subdivide it -- we can immediately get meshSpace's
   min and max vertex positions using the arrow operator, and the subdivision scale is set by the interaction range. 
6. Next, we'll define a configuration. The core thing the configuration stores is the arrangement of particles. You
   can set initial particle positions in one of two ways: 
	1. Creating a noise source and using setRandomParticlePositions, as is done in the sample code
	2. Feeding in your own positions. This can be done either as a set of meshPositions, if you already have 
	   barycentric coordinates, or, if you have R3 coordinates and are very confident they're on or nearly on
	   the mesh, you can use configuration->setMeshPositionsFromR3File. That takes a csv of rows, each with three entries,
	   that correspond the x, y, and z coordinates of the point in a traditional Cartesian representation.  
7. Having called setRandomParticlePositions(), we now define a force (of course, we could define
   the force, configuration, and updater in any order if we wanted). Here, we set a harmonicRepulsion force, 
   holding the energy scale (or stiffness) to unity and instead modifying it with our interaction range parameter, 
   which might be roughly thought of as particle size, via another make_shared shared pointer declaration. 
8. The force needs to be tied to a configuration so the force knows what to act on, so we use its setModel 
   function to tie the two together. 
9. We have a configuration and a force. In principle, we might define an updater next. In our example, though, 
   we instead create another shared pointer, this time to a simulation. 
10. Now that the simulation is defined, we can tie it together to the existing force function and configuration. 
   Remember that we could add any number of forces if we wanted, and they would all get evaluated at every 
   step. 
11. Finally, we define an updater. curvedSpaceSimulation.cpp allows for either NVE or gradientDescent, depending
    on the programBranch (just a way of specifying behavior). Let's assume our branch is either 0 (all-to-all 
    with no velocity) or 1 (cutoff radii and submeshed with no velocity). Then, we define our updater as a 
    gradientDescent object and tie it to the simulation via simulator->addUpdater(). 
12. Now, our simulator is tied to a configuration, force, and updater. With all of those defined, we create 
    our flat position vector (for efficient storage in a netcdf container), write a database object to store
    the results, and then call simulator->performTimestep() as many times as we want to advance the simulation 
    forward. This automatically calls the updater, which in turn calls the force and updates the configuration. 
  
   



