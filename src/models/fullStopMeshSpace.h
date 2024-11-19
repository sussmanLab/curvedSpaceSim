#ifndef fullStopMeshSpace_H
#define fullStopMeshSpace_H

#include "triangulatedMeshSpace.h"

class fullStopMeshSpace : public triangulatedMeshSpace
    {
    public:
        fullStopMeshSpace(){}; 

	//!Given a particle somewhere on the mesh, displace it in the direction of the vector, wrapping around faces; redefined from tms due to new boundary args
        virtual void transportParticleAndVectors(meshPosition &pos, vector3 &displacementVector, vector<vector3> &transportVectors);

    protected: 
	 //!Implement boundary conditions; these functions will force particles to stop if they hit a boundary of the space
        virtual void fullStopBoundaryEdge(pmpBarycentricCoordinates &targetBCs, pmpBarycentricCoordinates &sourceBCs, faceIndex &sourceFace, vertexIndex edgeV1, vertexIndex edgeV2, point3 innerVertex, vector<vector3> transportVectors, bool &continueShifting);
	virtual void fullStopBoundaryVertex(pmpBarycentricCoordinates &sourceBCs, pmpBarycentricCoordinates &targetBCs, vector3 &sourceNormal, faceIndex &sourceFace, vector3 &displacement, vector<point3> &vertexPositions, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vertexIndex intersectedV, bool &continueShifting);

	//!helper functions for boundary vertex behaviors -- these are designed to extend to more complex boundary conditions 
        virtual void projectVectorsForBoundaryVertex(vector3 heading, vertexIndex interesctedV, vertexIndex boundaryV, vector3 fNormal, faceIndex face, vector<vector3> &transportVectors);
        virtual void getBoundaryVertexHeading(vertexIndex v, faceIndex &sourceFace, vector3 &displacementVector, pmpBarycentricCoordinates &sourceBarycentricLocation, vector<point3> vertexPositions, vector3 sourceNormal, std::pair<vector3, vertexIndex> &result);

    };

#endif
