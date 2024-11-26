#ifndef tangentialMeshSpace_H
#define tangentialMeshSpace_H

#include "triangulatedMeshSpace.h"
#include "fullStopMeshSpace.h"

class tangentialMeshSpace : public fullStopMeshSpace 
    {
    public:
	//extends is the right word here -- we inherit fullstopmeshspace so we can use the helper functions for vertex intersections
        tangentialMeshSpace(){}; 

	//!Given a particle somewhere on the mesh, displace it in the direction of the vector, wrapping around faces; redefined from tms due to new boundary args
        virtual void transportParticleAndVectors(meshPosition &pos, vector3 &displacementVector, vector<vector3> &transportVectors);

    protected: 
	 //!Implement boundary conditions; these functions will direct particles along the boundary if they encounter one
        virtual void tangentialBoundaryEdge(pmpBarycentricCoordinates &sourceBCs, point3 &target, faceIndex &sourceFace, vector3 &sourceNormal, vertexIndex edgeV1, vertexIndex edgeV2, point3 innerVertex, vector<vector3> transportVectors, vector3 &displacement, halfedgeIndex &lastUsedHalfedge, bool &continueShifting);

	virtual void tangentialBoundaryVertex(pmpBarycentricCoordinates &sourceBCs, pmpBarycentricCoordinates &targetBCs, point3 &target, vector3 &sourceNormal, faceIndex &sourceFace, vector3 &displacement, vector<point3> &vertexPositions, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vertexIndex intersectedV, bool &continueShifting)
;

    };

#endif
