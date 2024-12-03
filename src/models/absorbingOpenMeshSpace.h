#ifndef absorbingOpenMeshSpace_H
#define absorbingOpenMeshSpace_H

#include "openMeshSpace.h"

class absorbingOpenMeshSpace : public openMeshSpace
    {
    public:
        absorbingOpenMeshSpace(){}; 

	//!Given a particle somewhere on the mesh, displace it in the direction of the vector, wrapping around faces; redefined from tms due to new boundary args
        virtual void transportParticleAndVectors(meshPosition &pos, vector3 &displacementVector, vector<vector3> &transportVectors);

    protected: 
	//!Implement boundary conditions; these functions will force particles to stop if they hit a boundary of the space
        virtual void boundaryEdge(pmpBarycentricCoordinates &targetBCs, pmpBarycentricCoordinates &sourceBCs, faceIndex &sourceFace, vertexIndex edgeV1, vertexIndex edgeV2, point3 innerVertex, vector<vector3> transportVectors, bool &continueShifting);
	virtual void boundaryVertex(pmpBarycentricCoordinates &sourceBCs, pmpBarycentricCoordinates &targetBCs, vector3 &sourceNormal, faceIndex &sourceFace, vector3 &displacement, vector<point3> &vertexPositions, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vertexIndex intersectedV, bool &continueShifting);
   
    };

#endif
