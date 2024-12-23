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
*/
class triangulatedMeshSpace : public baseSpace
    {
    public:
        triangulatedMeshSpace(){};

        //!load data from an off file and initialize all mesh data structures
        void loadMeshFromFile(std::string filename, bool verbose = false);

	//legacy declarations:
	/* 
        //!Given a particle somewhere on the mesh, displace it in the direction of the vector, wrapping around faces to make it a geodesic displacement
        virtual void displaceParticle(meshPosition &pos, vector3 &displacementVector, vector3 &force) = 0;
        //!On a mesh, if a particle goes across an edge the velocity vector gets rotated along with it
        virtual void transportParticleAndVelocity(meshPosition &pos, vector3 &v, vector3 &displacementVector, vector3 &force) = 0;
        */
	//!Given a particle somewhere on the mesh, displace it in the direction of the vector, wrapping around faces; any vectors given to the function are transported along with it
	virtual void transportParticleAndVectors(meshPosition &pos, vector3 &displacementVector, vector<vector3> &transportVectors); 


        //!Given a source particle and a vector of target points, determine the geodesic distance and store the start and end path tangents along the paths
        virtual void distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent, double distanceThreshold= VERYLARGEDOUBLE);
       
        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result);

	//!Get a random face and a random barycentric location within that face
        virtual void randomPosition(meshPosition &p, noiseSource &noise);

	//!Get a random face restricted to those included in a vector of faces and a random bary location within that face
	virtual void randomPositionWithinFaces(meshPosition &p, noiseSource &noise, vector<faceIndex> faces);

        //!Get a random unit vector in the tangent space of the targget position
        virtual void randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise);

	//!Functions to perform simple projections 'on' (align with) or 'off' (make perpendicular to) a given vector
	virtual void projectOn(vector3 &v, vector3 &direction);
	virtual void projectOff(vector3 &v, vector3 &direction);

	virtual void projectVectorsIfOverBoundary(vector<vector3> &vectors, vector3 orthogonal, vector3 inward);

	//!specialized function 
        //!Get the area of the mesh
        virtual double getArea(); 

        //!specialize the calculation of distances to use submeshes
        void distanceWithSubmeshing(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent,double distanceThreshold);

        //!Activate the use of submeshing routines
        void useSubmeshingRoutines(bool _useSubMesh, double maxDist = 1.0, bool _danger = false)
            {submeshingActivated = _useSubMesh;
            maximumDistance = maxDist;
            //dangerous submeshing no longer needed, but might be revived later TODO
            dangerousSubmeshing = _danger;
            };

        //!set a new cutoff scale for submeshing routines (used in some of Toler's mesh refinement tests)
    	void setNewSubmeshCutoff(double newCutoff)
            {
            useSubmeshingRoutines(true,newCutoff);
            }

	void printMaxDist() 
	    {
	    printf("Maximum distance is: %.6g \n", maximumDistance);
	    }

        void printMinMaxVP() 
            {
            printf("Min vertex position: (%.6g,%.6g,%.6g)  Max vertex position: (%.6g,%.6g,%.6g)  \n", minVertexPosition.x, minVertexPosition.y, minVertexPosition.z,  maxVertexPosition.x, maxVertexPosition.y, maxVertexPosition.z); 
	    }

        //!given a vector of meshPositions that represent barycentric coordinates, fill a second vector of meshPositions that represent the corresponding R3 positions
        void convertToEuclideanPositions(std::vector<meshPosition> &a, std::vector<meshPosition> &b);

	//!given a goal edge length, remesh the surface isotropically to have that edge length on average and set the new surface
        void isotropicallyRemeshSurface(double targetEdgeLength);


        //! The mesh object representing the surface
        triangleMesh surface;
        //! The lower-bottom-left position of the rectilinear domain containing the surface
        double3 minVertexPosition;
        //! The upper-top-right position of the rectilinear domain containing the surface
        double3 maxVertexPosition;        
	
	bool useTangentialBCs = false;


        void updateMeshSpanAndTree(bool updateSMSP = true);
    
    protected:
	void checkBaryNan(pmpBarycentricCoordinates bcoords, string message = "", int step = 0);
        void clampAndUpdatePosition(pmpBarycentricCoordinates &baryLoc, point3 &r3Loc, faceIndex &sFace, bool belowZero = false);
        //!Update some internal datastructures used in finding shortest paths
        //!Handle transport through vertices rather than across edges
        std::pair<faceIndex,vector3> throughVertex(vertexIndex &intersectedVertex, vector3 &toIntersection, faceIndex &sourceFace);
        bool verbose = false;
        shared_ptr<surfaceMeshShortestPath> globalSMSP;
        AABB_tree globalTree;
        bool submeshingActivated = false;
        //!A data structure for helping with submeshing routines
        submesher submeshAssistant;
        double maximumDistance = 0;
	//data structures associated with the potential for "dangerous" submeshes -- these can occur when submeshing routine is capable of returning a submesh with multiple connected components (e.g., when dealing with a surface that looks like an elephant's ear)
        bool dangerousSubmeshing = false;
        

        int maximumShiftEdgeCrossings = 1000000;
    };
#endif
