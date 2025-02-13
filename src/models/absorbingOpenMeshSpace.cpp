#include "absorbingOpenMeshSpace.h"
void absorbingOpenMeshSpace::updateShiftAtBoundaryEdge(
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
    bool& continueShifting)
    {
    if (transportVectors.size() > 0) 
        {
	//below handles just nonzero transport vectors size
        point3 ev1 = surface.point(edgeV1);
        point3 ev2 = surface.point(edgeV2);

        //check to see if any of our transported vectors were pointing over a boundary along with the displacement & project if so 
        vector3 orthogonalToEdge = normalize(CGAL::cross_product(PMP::compute_face_normal(sourceFace,surface), vector3(ev1,ev2)));
        vector3 inwardVector(globalSMSP->point(sourceFace, sourceBCs), innerVertex); 

	//if orthogonal is pointing out, then we want to project only if the overlap with the orthogonal vector is positive; 
        //if orthogonal is pointing in, then we want to project only if the overlap with the orthogonal vector is negative. 
        projectVectorsIfOverBoundary(transportVectors, orthogonalToEdge, inwardVector);
	}
    //this is all that the boundary conditions actually do -- stop the motion and move the particle to the edge that it intersected.     
    continueShifting=false; 
    targetBCs = sourceBCs;  
    }

void absorbingOpenMeshSpace::updateShiftAtBoundaryVertex(
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
    bool& continueShifting)
    {
    //first get heading & next vertex along the boundary in that heading; source face is modified
    std::pair<vector3, vertexIndex> headingAndVertex;
    getBoundaryVertexHeading(intersectedV, sourceFace, displacement, sourceBCs, vertexPositions, sourceNormal, headingAndVertex);
    //project any vectors we need to project according to the calculated boundary heading
    if (transportVectors.size() > 0)
        projectVectorsForBoundaryVertex(headingAndVertex.first, intersectedV, headingAndVertex.second, sourceNormal, sourceFace, transportVectors);
    
    //technically, the above moves us into an adjacent face -- since the vertex is in the same place
    //either way, it doesn't change the outcome, but this version is more extensible to tangential
    //BCs at the cost of more code overhead.  
    continueShifting = false; 
    targetBCs = sourceBCs; //soruce bcs are already intersection point -- this guarantees we don't accidentally move after loop is done
     
    halfedgeIndex nullHalfedge(-1); 
    lastUsedHalfedge = nullHalfedge;  
    }




