#ifndef triangulatedMeshSpace_H
#define triangulatedMeshSpace_H

#include "functionUtilities.h"
#include "meshUtilities.h"
#include "baseSpace.h"
#include "submesher.h"
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

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
        //!On a mesh, if a particle goes across an edge the velocity vector gets rotated along with it
        virtual void transportParticleAndVelocity(meshPosition &pos, vector3 &v, vector3 &displacementVector);

        //!Given a source particle and a vector of target points, determine the geodesic distance and store the start and end path tangents along the paths
        virtual void distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent, double distanceThreshold= VERYLARGEDOUBLE);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result);

        virtual void randomPosition(meshPosition &p, noiseSource &noise);

        virtual void randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise);

        void distanceWithSubmeshing(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent,double distanceThreshold);

        void useSubmeshingRoutines(bool _useSubMesh, double maxDist = 1.0, bool _danger = false)
            {submeshingActivated = _useSubMesh;
            maximumDistance = maxDist;
            //dangerous submeshing no longer needed, but might be revived later
            dangerousSubmeshing = _danger;
            };

	void setNewSubmeshCutoff(double newCutoff)
        {
           if(submeshingActivated) maximumDistance = newCutoff;
           else
           {
               std::cout << "submeshing is not activated, cannot set new cutoff" <<std::endl;
               throw std::exception();
           }
        }

        //!given a vector of meshPositions that represent barycentric coordinates, fill a second vector of meshPositions that represent the corresponding R3 positions
        void convertToEuclideanPositions(std::vector<meshPosition> &a, std::vector<meshPosition> &b);

	//!given a goal edge length, remesh the surface isotropically to have that edge length on average and set the new surface
        void isotropicallyRemeshSurface(double targetEdgeLength);

        //data structures
        triangleMesh surface;
        double3 minVertexPosition;
        double3 maxVertexPosition;
    protected:
        void updateMeshSpanAndTree();
        std::pair<faceIndex,vector3> throughVertex(vertexIndex &intersectedVertex, vector3 &toIntersection, faceIndex &sourceFace);
	bool verbose = false;
        shared_ptr<surfaceMeshShortestPath> globalSMSP;
        AABB_tree globalTree;
        //data structures associated with submeshing routines
        bool submeshingActivated = false;
        submesher submeshAssistant;
        double maximumDistance = 0;
        //data structures associated with the potential for "dangerous" submeshes -- these can occur when submeshing routine is capable of returning a submesh with multiple connected components (e.g., when dealing with a surface that looks like an elephant's ear)
        bool dangerousSubmeshing = false;

        int maximumShiftEdgeCrossings = 1000000;
    };
#endif
