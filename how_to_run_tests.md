#Test Explainer
Each of the tests is essentially run through command line inputs, 
parsed via TCLAP. They're all variations on the core sample code 
in curvedSpaceSimulation.cpp. For each of the time vs. particle
number and time vs. mesh complexity tests, we strictly 
utilized harmonically repulsive forces for flexible numbers of 
iterations at each value tested.

Each of the below should be executable once you have successfully
invoked 
```
cmake ..
make 
``` 
in a build directory. 

##Cost vs. N 
As with all of these, one should be able to run the cost vs. 
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

The command line input to remake the trend in our figure 3 
using a reproducible rng would be: 
```
./t_v_n.out -i 5 -n 15 -m "../exampleMeshes/torus_isotropic_remesh.off" -a 1 -t .01 -r
```

This input will run both all-to-all and submeshed tests
using the interaction range specified by -a. The all-to-all
tests take quite a long time (and in fact, they are the limiting
factor on the particle number scale we can use for these tests -- 15 
was about as many doublings as our testing machine could take and
finish in reasonable time), so if one only wishes to reproduce the scaling behavior
for submeshed tests, please comment out the all-to-all test
at the designated comment markers. One file will be produced
with the particle numbers and paired average simulation timestep 
costs. 

To reproduce the trends in the inset, one needs only specify the interaction range as 
.5, .25, and .125 in successsion and plot those results alongside the 
innteraction range 1 result.  

##Cost vs. Complexity

The complexity tests rely on successive remeshings of the original mesh. 
To do this, we utilize CGAL's Polygon Mesh Processing's isotropic 
remesher, which remeshes toward a target mean edge length after remeshing. 
The remesher is built into our source code explicitly for this test, 
but can behave somewhat strangely when placed under stress; namely,
we observed some strange indexing behavior that could, in very specific
scenarios, cause the intersection check of shift to fail
(namely when a source and target were both within numerical precision of an edge).  

For our test, we remesh the surface, then write the remeshed 
surface to a file, then read it. This routine completely avoids the
indexing behavior that caused the above error, but is a little slower. 

If one has already generated remeshes, which we provide in the remeshes folder,
then these can be loaded in to perform the test. E.g., to replicate the submeshed
elephant scaling, one can run the command:  

```
./t_v_complexity.out -m "../exampleMeshes/triangulatedElephant.off" -w "elephant" -z 1 -y 0
```
where similar commands with the writeFile & base mesh changed can reproduce the 
sphere, 3-1 torus, and 20-1 torus results: 
```
./t_v_complexity.out -m "../exampleMeshes/sphere_rb20.off" -w "sphere" -z 1 -y 0
./t_v_complexity.out -m "../exampleMeshes/big_torus_rb50.off" -w "big_torus" -z 1 -y 0
./t_v_complexity.out -m "../exampleMeshes/torusrb20.off" -w "torus" -z 1 -y 0
```

Here, -y 0 indicates that the code is not to write new remeshes, and instead use 
the existing ones in the remeshes folder. The submeshes that 
we used in our tests have a lower bound mean edge length of 1/4 of the original; 
that's because extremely complex meshes require extremely complex sequence trees,
and those huge sequence trees could cause the all-to-all tests to fail. Of course,
the submeshed tests should be resilient to much more complex meshes. If one wants
to investigate much higher complexity remeshes, the line to edit is line 155, 
near the end of the cpp file and marked with a comment.  

##Torus Relaxations

Because the torus relaxations are the closest to a true numerical physics
experiment of the tests, the torusCrystallization cpp file has the largest set of command 
line options. One can set  the number of particles to be relaxed, the force
tolerance at which to stop crystallization, the harmonic stiffness, and a 
file for starting locations in addition to the usual choice of mesh, 
interaction range, reproducibility, submeshing on/off, 
and save frequency. We always had the reproducibility switch set to default (on) 
for our tests. To see which command line option refers to which thing, 
please refer to the torusCrystallization.cpp file.  

The command we used to produce the 32 particle torus configuration
shown in the paper was 
```
./torusCrystallization.out -m "../exampleMeshes/torusrb20.off" -n 32 -s 100 -i 7 -e .01 -z 1 -a 2.25
```

For 1000 particles, we used the command
```
./torusCrystallization.out -m "../exampleMeshes/torusrb20.off" -n 1000 -s 100 -i 7 -e .01 -z 1 -a 0.4
``` 
In both of the above, the final configuration is very sensitive to the interaction radius.
E.g., one should expect very different results for -a .41 in the 1000 particle test 
vs. -a .4. We performed tests not included in the paper at, for example -a .43 that still produced
realistic results. 
For each, we chose interaction radii that corresponded to disks that would approximately
cover at or above 100% of the torus' surface area -- however, true toroidal packing with all disks 
only kissing would instead correspond to around 90% coverage. 

Reproducing the Voronoi diagrams is much harder, and requires, at present, a Windows machine. 
We used the SurfaceVoronoi package, located at https://github.com/sssomeone/SurfaceVoronoi.
To generate bisectors, we ran the above within Windows VisualStudio using the torusrb20.off
mesh as the surface and a 
```
[particle number]_[fn cutoff exponent]_particle_voronoi_positions.csv
```
file as position input. These files are combinations of face indices and r3 positions for each particle 
at the end of relaxation. We specifically used the SurfaceVoronoi function GetLRVDBisectorsWithDui
to generate both a bisector image file and a triangular connectivity object file. 
We will be happy to provide our specific edited SurfaceVoronoi installation upon 
request, although it is not part of this codebase. 


