#ifndef openMeshSpace_H
#define openMeshSpace_H

#include "triangulatedMeshSpace.h"

class openMeshSpace : public triangulatedMeshSpace
    {
    public: 
        openMeshSpace(){};

	virtual void transportParticleAndVectors(meshPosition &pos, vector3 &displacementVector, vector<vector3> &transportVectors);
     //! functions signature for implementation by child classes. If this signature does not include an argument necessary for your implementation, add the relevant argument. 
        virtual void updateAtBoundaryEdge(
    pmpBarycentricCoordinates& sourceBCs,
    pmpBarycentricCoordinates& targetBCs,
    point3& target,
    faceIndex& sourceFace,
    vector3& sourceNormal,
    const vertexIndex edgeV1,
    const vertexIndex edgeV2,
    const point3& innerVertex,
    vector<vector3>& transportVectors,
    vector3& displacement,
    halfedgeIndex& lastUsedHalfedge,
    bool& continueShifting
)     
        {
        ERRORERROR("Boundary edge: base open mesh space does not have boundary conditions."); 
        }

virtual void updateAtBoundaryVertex(
    pmpBarycentricCoordinates& sourceBCs,
    pmpBarycentricCoordinates& targetBCs,
    point3& target,
    vector3& sourceNormal,
    faceIndex& sourceFace,
    vector3& displacement,
    const vertexIndex intersectedV,
    const vector<point3>& vertexPositions,
    vector<vector3>& transportVectors,
    halfedgeIndex& lastUsedHalfedge,
    bool& continueShifting
)
        {
	ERRORERROR("Boundary vertex: base open mesh space does not have boundary conditions."); 
        }
    
    protected:
	//!helper functions for boundary vertex behaviors -- these are designed to extend to more complex boundary conditions
        virtual void projectVectorsForBoundaryVertex(vector3 heading, vertexIndex interesctedV, vertexIndex boundaryV, vector3 fNormal, faceIndex face, vector<vector3> &transportVectors);
        virtual void getBoundaryVertexHeading(vertexIndex v, faceIndex &sourceFace, vector3 &displacementVector, pmpBarycentricCoordinates &sourceBarycentricLocation, vector<point3> vertexPositions, vector3 sourceNormal, std::pair<vector3, vertexIndex> &result);

    };

#endif


