#ifndef triangulatedMeshSpace_H
#define triangulatedMeshSpace_H

#include "functionUtilities.h"
#include "meshUtilities.h"
#include "baseSpace.h"

/*! \file triangulatedMeshSpace.h"
* \brief defines an interface to CGAL mesh-based functionality
*/

//! A class that interfaces with CGAL functionality

/*!
On "loadMeshFromFile", loads a triangle mesh into the surface member and also initializes a global
surfaceMeshShortestPath (and builds an associated AABB_tree), which can be used for any whole-mesh
operations
The triangulatedMeshSpace *ASSUMES* that the data contained in the meshPositions is not actually a
(point3,faceIndex) pair, but rather a point3 which is storing the 3 Barycentric coordinates of the
point. (and a faceIndex integer)
*/
class triangulatedMeshSpace : public baseSpace
    {
    public:
        triangulatedMeshSpace(){};

        //!load data from an off file and initialize all mesh data structures
        void loadMeshFromFile(std::string filename, bool verbose = false);


        //!Given a particle somewhere on the mesh, displace it in the direction of the vector, wrapping around faces to make it a geodesic displacement
        virtual void displaceParticle(meshPosition &pos, vector3 &displacementVector);

        //!Given a source particle and a vector of target points, determine the geodesic distance and store the start and end path tangents along the paths
        virtual void distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result);

        void distanceWithSubmeshing(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent);

        void useSubmeshingRoutines(bool _useSubMesh){submeshingActivated = _useSubMesh;};

        //!given a vector of meshPositions that represent barycentric coordinates, fill a second vector of meshPositions that represent the corresponding R3 positions
        void convertToEuclideanPositions(std::vector<meshPosition> &a, std::vector<meshPosition> &b);

        //data structures
        triangleMesh surface;
    protected:
        bool submeshingActivated = false;
        shared_ptr<surfaceMeshShortestPath> globalSMSP;
        AABB_tree globalTree;
    };
#endif
