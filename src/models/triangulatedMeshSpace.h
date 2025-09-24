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

//! A class that interfaces with CGAL functionality for working with triangulated meshes
/*!
Note that this class interfaces tightly with some of the function in the meshUtilities files
Currently, the code always considers boundary conditions -- for closed surfaces, they just aren't called. If useTangentialBCs is on, 
transportParticleAndVectors will slide particles along the boundary rather than stop them. Future versions of the code will differentiate
routines that use boundary conditions at all with routines that don't (for solely closed surfaces) in a modular way. 
*/
class triangulatedMeshSpace : public baseSpace
    {
    public:
        triangulatedMeshSpace(){};

        //!load data from an off file and initialize all mesh data structures
        virtual void loadMeshFromFile(std::string filename, bool verbose = false);

        //!displace particle just calls the transportParticleAndVectors function
        virtual void displaceParticle(meshPosition &pos, vector3 &displacementVector);
        
        //!Given a particle somewhere on the mesh, displace it in the direction of the vector, wrapping around faces; any vectors given to the function are transported along with it
        virtual void transportParticleAndVectors(meshPosition &pos, vector3 &displacementVector, vector<vector3> &transportVectors);

        //!Given a source particle and a vector of target points, determine the geodesic distance and store the start and end path tangents along the paths
        virtual void distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent, double distanceThreshold= VERYLARGEDOUBLE);

        //!Get a random face and a random barycentric location within that face
        virtual void randomPosition(meshPosition &p, noiseSource &noise);

        //!Get a random unit vector in the tangent space of the targget position
        virtual void randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise);

        //!Get the area of the mesh
        virtual double getArea();

        //!specialize the calculation of distances to use submeshes
        void distanceWithSubmeshing(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent,double distanceThreshold);

        //!Activate the use of submeshing routines
        void useSubmeshingRoutines(bool _useSubMesh, double maxDist = 1.0, bool _danger = false)
            {
            submeshingActivated = _useSubMesh;
            maximumDistance = maxDist;
            //dangerous submeshing no longer needed, but might be revived later TODO
            dangerousSubmeshing = _danger;
            };

        //!set a new cutoff scale for submeshing routines (used in some of Toler's mesh refinement tests)
    	void setNewSubmeshCutoff(double newCutoff)
           {
           useSubmeshingRoutines(true,newCutoff);
           }

        //!given a vector of meshPositions that represent barycentric coordinates, fill a second vector of meshPositions that represent the corresponding R3 positions
        void convertToEuclideanPositions(std::vector<meshPosition> &a, std::vector<meshPosition> &b);

        //!given a goal edge length, remesh the surface isotropically to have that edge length on average and set the new surface
        void isotropicallyRemeshSurface(double targetEdgeLength);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result);

        //!The actual meshed surface itself
        triangleMesh surface;
        //!The lower-bottom-left position of the rectilinear domain containing the surface
        double3 minVertexPosition;
        //!The upper-top-right position of the rectilinear domain containing the surface
        double3 maxVertexPosition;
        
	virtual void printSourceTargetDisplacementInfo(point3 sourcePoint, point3 target, pmpBarycentricCoordinates sourceBarycentricLocation, pmpBarycentricCoordinates targetBarycentricLocation, vector3 displacementVector);

    protected:
        //!Update some internal datastructures used in finding shortest paths
        void updateMeshSpanAndTree();
        //!Handle transport through vertices rather than across edges
        std::pair<faceIndex,vector3> throughVertex(vertexIndex &intersectedVertex, vector3 &toIntersection, faceIndex &sourceFace);
        //!Helper functions to handle different triangle border crossing cases
        virtual void updateForEdgeIntersection(pmpBarycentricCoordinates &sourceBarycentricLocation, point3 &sourcePoint, pmpBarycentricCoordinates &intersectionPoint, vector3 &currentSourceNormal, faceIndex &currentSourceFace, point3 &target, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vector3 &displacementVector, vector<point3> &vertexPositions, halfedgeIndex intersectedEdge);

        virtual void updateForVertexIntersection(pmpBarycentricCoordinates &sourceBCs, point3 &sourcePoint, faceIndex &sourceFace, point3 &target, vector3 &displacementVector, vector3 &sourceNormal, vector<point3> vertexPositions, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vertexIndex intersectedV, vector3 toIntersection);

        bool verbose = false;
        shared_ptr<surfaceMeshShortestPath> globalSMSP;
        AABB_tree globalTree;
        bool submeshingActivated = false;
        //!specialized function to assist with boundary conditions -- this is an outcome subfunction for boundary behaviors
        virtual void projectVectorsIfOverBoundary(vector<vector3> &vectors, vector3 orthogonal, vector3 inward);
        //!A data structure for helping with submeshing routines
        submesher submeshAssistant;
        double maximumDistance = 0;
        //data structures associated with the potential for "dangerous" submeshes -- these can occur when submeshing routine is capable of returning a submesh with multiple connected components (e.g., when dealing with a surface that looks like an elephant's ear)
        bool dangerousSubmeshing = false;
        int maximumShiftEdgeCrossings = 1000000;
    };
#endif
