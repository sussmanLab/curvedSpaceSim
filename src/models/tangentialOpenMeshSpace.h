#ifndef tangentialOpenMeshSpace_H
#define tangentialOpenMeshSpace_H

#include "openMeshSpace.h"

class tangentialOpenMeshSpace : public openMeshSpace 
    {
    public:
        tangentialOpenMeshSpace(){}; 

    protected: 
	 //!Implement boundary conditions; these functions will direct particles along the boundary if they encounter one
         virtual void updateShiftAtBoundaryEdge(pmpBarycentricCoordinates& sourceBCs, 
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
    bool& continueShifting);
         
	 virtual void updateShiftAtBoundaryVertex(
    pmpBarycentricCoordinates& sourceBCs, 
    pmpBarycentricCoordinates& targetBCs, 
    point3& target, 
    vector3& sourceNormal, 
    faceIndex& sourceFace, 
    vector3& displacement, 
    vertexIndex intersectedV, 
    vector<point3>& vertexPositions, 
    vector<vector3>& transportVectors, 
    halfedgeIndex& lastUsedHalfedge, 
    bool& continueShifting);
    };

#endif
