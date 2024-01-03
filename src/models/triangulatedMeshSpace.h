#ifndef triangulatedMeshSpace_H
#define triangulatedMeshSpace_H

#include "baseSpace.h"

/*! \file triangulatedMeshSpace.h"
* \brief defines an interface to CGAL mesh-based functionality
*/

//! A class that interfaces with CGAL functionality


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

        void distanceWithSubmeshing(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent);

        void useSubmeshingRoutines(bool _useSubMesh){submeshingActivated = _useSubMesh;};

        void findIntersection(faceIndex sourceFace, point3 source, point3 target, std::vector<vertexIndex> &vertexIndices, point3 &intersectionPoint, std::vector<vertexIndex> &intersections);

        //data structures
        triangleMesh surface;
    protected:
        bool submeshingActivated = false;

        double vectorMagnitude(vector3 v);
    };
#endif
