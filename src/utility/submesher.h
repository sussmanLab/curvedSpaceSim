#ifndef submesher_H
#define submesher_H

#include "meshUtilities.h"

/*! \file submesher.h */

/*!
A class that assists with sub-mesh operation... given a triangleMesh and a set of
target points, return a local patch of  the mesh which contains all of the points
and the intermediate faces needed to determine local geodesic paths.
Currently I am making this a class in case I want to include accelerating data 
structures later.
*/
class submesher
    {
    public:
        submesher(){};

        triangleMesh constructSubmeshFromSourceAndTargets(triangleMesh &mesh, faceLocation &source, std::vector<faceLocation> &targets, double &maximumDistanceFromSource);
    };
#endif
