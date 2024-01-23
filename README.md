# curvedSpaceSimulation

curvedSpaceSimulation implements molecular-dynamics-like evolution of degrees of
freedom on meshed (triangulated) surfaces embedded in 3D space. All motion of
particles follows discrete geodesics on a given surface, and all distances are compute
as geodesic distances along the surface (i.e., particles interact via "curved lines
of force" rather than interacting according to the Euclidean distance). Parallel transport (e.g., of velocity vectors during the update of equations of motion) is 
implemented as part of what it means to displace a particle along a geodesic.

The standard all-to-all geodesic-finding algorithm implemented has O(v^2) computational
complexity, where v is the number of vertices of the surface mesh (although, in
practice, the scaling does not seem to be that bad on representative surfaces). In the
case of simulations of particles with short-range interactions, submeshing routines
have been implemented so that the computational cost scales like N log N (where N is
the particle number) and like (w^2), where w is the typical number of vertices in a patch
of surface spanning the maximum interaction range.

Information on installing the project can be found [here](/Installation.md). A very 
rough outline of some of the main classes andthe basic operating flow of the primary
branches can be found [here](/doc/BasicInformation.md). A very brief description of
the main executables distributed can be found [here](/doc/sampleCode.md).

## Project documentation

Additional files in the doc director provide some basic documentation

[Basic class overview](/doc/BasicInformation.md)

[Installation guide](/Installation.md)

[Sample code snippets](/doc/sampleCode.md)

[Citations](/Citations.md)

[Contributors](/Contributors.md)

[Change log](/Changelog.md)

[Open-source information](/LICENSE)
