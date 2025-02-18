#include "tangentialOpenMeshSpace.h"

void tangentialOpenMeshSpace::updateAtBoundaryEdge(
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
    point3 ev1 = surface.point(edgeV1);
    point3 ev2 = surface.point(edgeV2);
    vector3 edgeVectorForward(ev1,ev2); 
    vector3 edgeVectorBackward(ev2,ev1);
    edgeVectorForward = normalize(edgeVectorForward); 
    edgeVectorBackward = normalize(edgeVectorBackward); 
    double forwardDot = normalize(displacement)*edgeVectorForward;
    double backwardDot = normalize(displacement)*edgeVectorBackward;
    double displacementLength = vectorMagnitude(displacement); 

    //before we update displacement, check if any transported vectors point over a boundary 
    vector3 orthogonalToEdge = normalize(CGAL::cross_product(sourceNormal, vector3(ev1,ev2))); 
  
    //if orthogonal is pointing out, project only if the overlap with it is positive
    //if orthogonal is pointing in, project only if the overlap with it is negative
    vector3 inwardVector(globalSMSP->point(sourceFace, sourceBCs), innerVertex); 
    projectVectorsIfOverBoundary(transportVectors, orthogonalToEdge, inwardVector);

    //now that the displacement vector has determined whether to project, 
    //update displacement vector 
    if ((forwardDot <= 0) && (backwardDot <= 0))
        continueShifting = false;
    if (forwardDot > backwardDot) 
        displacement = forwardDot*displacementLength*edgeVectorForward;
    else 
        displacement = backwardDot*displacementLength*edgeVectorBackward; 
    
    target = globalSMSP->point(sourceFace, sourceBCs) + displacement; 
    halfedgeIndex nullHalfedge(-1);
    lastUsedHalfedge = nullHalfedge; 
    }

void tangentialOpenMeshSpace::updateAtBoundaryVertex(
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
    
    //getBoundaryVertexHeading updates the source face, barycentric coordinates, vertex positions, and normal, 
    //so all we need to do is point our vector along the boundary heading that was output
    projectVectorOntoDirection(displacement, headingAndVertex.first);
    //source is already updated to be at the intersection point
    target = globalSMSP->point(sourceFace, sourceBCs)+displacement;
    
    halfedgeIndex nullHalfedge(-1); 
    lastUsedHalfedge = nullHalfedge;  
    }
