#ifndef submesher_H
#define submesher_H

#include "meshUtilities.h"
#include <set>
typedef PMP::Face_location<triangleMesh, FT>                            faceLocation;

/*! \file submesher.h */

/*!
A class that assists with sub-mesh operation... given a triangleMesh and a set of
target points, return a local patch of the mesh which contains all of the points
and the intermediate faces needed to determine local geodesic paths.
Currently I am making this a class in case I want to include accelerating data 
structures later.
*/
class submesher
    {
    public:
        submesher(){};

        //Note that the faceMap gets filled with a map to translate between the original faceIndex of the mesh and the new faceIndex of the submesh
        //! create a triangle mesh and return it. It will be a subset of the original mesh that contains the vertices and faces of the source and all target particles within some maximum distance. 
        triangleMesh constructSubmeshFromSourceAndTargets(triangleMesh &mesh, faceLocation &source, std::vector<faceLocation> &targets, double &maximumDistanceFromSource, std::unordered_map<vertexIndex,int> &vertexMap, std::unordered_map<faceIndex,int> &faceMap);

    protected:
        //!A helper function for building a submesh given a collection of faces
        triangleMesh constructSubmeshFromFaceSet(triangleMesh &mesh, std::unordered_set<faceIndex> &faces,std::unordered_map<vertexIndex,int> &vertexMap, std::unordered_map<faceIndex,int> &faceMap);
    };
#endif
