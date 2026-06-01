#ifndef reflectingOpenMeshSpace_H
#define reflectingOpenMeshSpace_H

#include "openMeshSpace.h"

/*!
reflectingOpenMeshSpace implements boundary conditions relevant when a mesh has an open
boundary. When a particles hits a boundary, all transported vectors have their components
normal to the boundary reversed, and in that timestep the particle will travel away from
the boundary at the speed it collided with it. 
*/

class reflectingOpenMeshSpace : public openMeshSpace
{
public:
    reflectingOpenMeshSpace(){};

protected:
    //! Implement boundary conditions; these functions will direct particles along the boundary
    //! if they encounter one & redirect any transported vectors to remove components pointing over the boundary
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
        bool& continueShifting);

    virtual void updateAtBoundaryVertex(
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

    //! Reflect transported vectors at a boundary vertex
    virtual void reflectVectorsForBoundaryVertex(
        vector3 heading,
        vertexIndex intersectedV,
        vertexIndex boundaryV,
        vector3 fNormal,
        faceIndex face,
        vector<vector3>& transportVectors);

    //! Modify any vectors that cross a boundary by reversing the boundary-orthogonal component
    virtual void reflectVectorsIfOverBoundary(
        vector<vector3>& vectors,
        vector3 orthogonal,
        vector3 inward);
};

#endif
