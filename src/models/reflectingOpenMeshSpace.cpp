#include "reflectingOpenMeshSpace.h"

void reflectingOpenMeshSpace::reflectVectorsIfOverBoundary(vector<vector3> &vectors, vector3 orthogonal, vector3 inward)
    {
    double inOutDot = orthogonal*inward;
    bool pointsOut = (inOutDot < 0);
    int nVectors = vectors.size();

    for (int i = 0; i < nVectors; i++)
        {
        vector3 v = vectors[i];
        double vAndPerpDot = v*orthogonal;
        bool vAlongPerp = (vAndPerpDot > 0);
        //if both the velocity is along the perpendicular direction AND the perpendicular director points out
	// or BOTH are not true, then reflect the vector
	if ((pointsOut && vAlongPerp) || (!pointsOut && !vAlongPerp))
	    {
	    vector3 vPerp = (normalize(orthogonal)*v)*normalize(orthogonal);
            v = v - 2*vPerp;
	    }
	vectors[i] = v;
        }
    }

void reflectingOpenMeshSpace::reflectVectorsForBoundaryVertex(vector3 heading, vertexIndex intersectedV, vertexIndex boundaryV, vector3 fNormal, faceIndex face, vector<vector3> &transportVectors)
    {
    //now project vectors iff they are pointing over the new boundary           
    vector3 orthogonalToEdge = normalize(CGAL::cross_product(fNormal, heading));
    vector<vertexIndex> faceVertices;
    faceVertices.push_back(vertexIndex(0));
    faceVertices.push_back(vertexIndex(0));
    faceVertices.push_back(vertexIndex(0));

    getVertexIndicesFromFace(surface,face, faceVertices);

    point3 insideVertex(0,0,0);

    for (int i = 0; i < 3; i ++)
        {
        vertexIndex sourceV = faceVertices[i];
        if ((sourceV == intersectedV) || (sourceV == boundaryV))
            continue;
        insideVertex = surface.point(sourceV);
        }
        vector3 inwardVector(surface.point(intersectedV), insideVertex);
        reflectVectorsIfOverBoundary(transportVectors, orthogonalToEdge, inwardVector);
    }


void reflectingOpenMeshSpace::updateAtBoundaryEdge(
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
    //catalog the vertices of the edge we intersected & which direction we're going via dots
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
    reflectVectorsIfOverBoundary(transportVectors, orthogonalToEdge, inwardVector);

    //now that the displacement vector has determined whether to project, 
    //update displacement vector 
    vector3 projection(0,0,0);
    vector3 perpendicular(0,0,0);
    if ((forwardDot <= 0) && (backwardDot <= 0))
        continueShifting = false;
    if (forwardDot > backwardDot) {
        projection = forwardDot*edgeVectorForward;
        perpendicular = displacement - projection; 
	displacement = displacementLength*normalize(projection - perpendicular);
    	}
    else { 
        projection = backwardDot*edgeVectorBackward;
        perpendicular = displacement - projection;
        displacement = displacementLength*normalize(projection - perpendicular);	
    	}
    
    target = globalSMSP->point(sourceFace, sourceBCs) + displacement; 
    halfedgeIndex nullHalfedge(-1);
    lastUsedHalfedge = nullHalfedge; 
    }

void reflectingOpenMeshSpace::updateAtBoundaryVertex(
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
        reflectVectorsForBoundaryVertex(headingAndVertex.first, intersectedV, headingAndVertex.second, sourceNormal, sourceFace, transportVectors);
    

    //getBoundaryVertexHeading updates the source face, barycentric coordinates, vertex positions, and normal.
    //This function uses the heading found as the basis for a reflection. We have to pick either the source edge
    //or next edge; here the next edge is chosen since we'll be traveling in the next face.
    vector3 projection = (displacement*headingAndVertex.first)*headingAndVertex.first;
    vector3 perpendicular = displacement - projection; 
    displacement = vectorMagnitude(displacement)*normalize(projection - perpendicular);

    //source is already updated to be at the intersection point
    target = globalSMSP->point(sourceFace, sourceBCs)+displacement;
    
    halfedgeIndex nullHalfedge(-1); 
    lastUsedHalfedge = nullHalfedge;  
    }
