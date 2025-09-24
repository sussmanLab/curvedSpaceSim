#ifndef absorbingOpenMeshSpace_H
#define absorbingOpenMeshSpace_H

#include "openMeshSpace.h"
/*!
absorbingOpenMeshSpace implements boundary conditions relevant when a mesh has an open
boundary. When a particles hits a boundary a particle stops where it hits the edge
and any transported vectors (if any) have their component normal to the boundary set to zero
*/
class absorbingOpenMeshSpace : public openMeshSpace
    {
    public:
        absorbingOpenMeshSpace(){}; 

    protected: 
        //!Implement boundary conditions; these functions will stop a particle if it hits a boundary. continueShifting is set to false, vectors are transported to be in-surface, target point is set to boundary intersection. 
        void updateAtBoundaryEdge(pmpBarycentricCoordinates& sourceBCs, 
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
         
        void updateAtBoundaryVertex(
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
