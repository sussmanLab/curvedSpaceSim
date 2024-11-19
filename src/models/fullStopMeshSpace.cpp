#include "fullStopMeshSpace.h"

void fullStopMeshSpace::getBoundaryVertexHeading(vertexIndex v, faceIndex &sourceFace, vector3 &displacementVector, pmpBarycentricCoordinates &sourceBarycentricLocation, vector<point3> vertexPositions, vector3 sourceNormal, std::pair<vector3, vertexIndex> &result) 
    {
    std::vector<faceIndex> faceIndices;
    std::vector<std::pair<vertexIndex,faceIndex>> neighborsAndFaces; 
    faceIndices.reserve(8);
    neighborsAndFaces.reserve(8); 

    faceIndex cand;	
    faceCirculator fbegin(surface.halfedge(v), surface), fdone(fbegin);
    do 
        {
        cand = *fbegin++; 
        if (cand == sourceFace) 
            continue;
        else 
            faceIndices.push_back(cand); 
        } while (fbegin != fdone);       

    vertexCirculator vbegin(surface.halfedge(v),surface), done(vbegin);
    do 
        {
        for (faceIndex f: faceIndices) 
            { 
            if (surface.face(surface.halfedge(*vbegin++)) == f)
                neighborsAndFaces.push_back(make_pair(*vbegin++, f)); 
            }	
        } while(vbegin != done);	

    //below serves the rule that the next heading will be along the edge that maximally 
    //overlaps with the present heading. 
    double maxOverlap = -1; //dummy value, should always be overwritten
    double currentOverlap = 0; 
    vector3 boundaryHeading = vector3(0,0,0); 
    faceIndex nextFace = sourceFace;
    vertexIndex boundaryVertex = v; 

    for (auto vi: neighborsAndFaces) 
        {
        vector3 outwardVector(surface.point(v),surface.point(vi.first));
        outwardVector = normalize(outwardVector);
        currentOverlap = outwardVector*normalize(displacementVector); 
        if (currentOverlap > maxOverlap) 
            {
            boundaryHeading = outwardVector; 
            maxOverlap = currentOverlap;
            nextFace = vi.second;
            boundaryVertex = vi.first;
            }
        }

    if (nextFace == sourceFace)
        ERRORERROR("New face not found for boundary vertex crossing.");
    if (boundaryVertex == v) 
        ERRORERROR("vertex not updated in throughvertex boundary crossing");
   
    // last step is doing small updates to source face, which happen by address here 
    // as it is the most convenient place to do so -- technically none of this is necessary
    // for full stop  
    sourceFace = nextFace;
    getVertexPositionsFromFace(surface, sourceFace, vertexPositions);
    sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0], vertexPositions[1], vertexPositions[2], surface.point(v));
    sourceNormal = PMP::compute_face_normal(sourceFace, surface);   
    result.first = boundaryHeading;
    result.second = boundaryVertex; 
    }

void fullStopMeshSpace::projectVectorsForBoundaryVertex(vector3 heading, vertexIndex intersectedV, vertexIndex boundaryV, vector3 fNormal, faceIndex face, vector<vector3> &transportVectors) 
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
        projectVectorsIfOverBoundary(transportVectors, orthogonalToEdge, inwardVector);	
    }


void fullStopMeshSpace::fullStopBoundaryVertex(pmpBarycentricCoordinates &sourceBCs, pmpBarycentricCoordinates &targetBCs, vector3 &sourceNormal, faceIndex &sourceFace, vector3 &displacement, vector<point3> &vertexPositions, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vertexIndex intersectedV, bool &continueShifting)
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

void fullStopMeshSpace::fullStopBoundaryEdge(pmpBarycentricCoordinates &targetBCs, pmpBarycentricCoordinates &sourceBCs, faceIndex &sourceFace, vertexIndex edgeV1, vertexIndex edgeV2, point3 innerVertex, vector<vector3> transportVectors, bool &continueShifting)
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


/*
Core displacement routine. Takes a degree of freedom and a direction + distance to move in, and moves 
it around the mesh, always remaining in the tangent plane of a mesh face. Any vectors passed in as transportVectors 
are parallel transported along the path of the particle.

As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex). 

pos: a meshPosition object (three doubles representing barycentric coordinates and a face index) 
     representing the position of a degree of freedom. 
displacementVector: The (euclidean) direction and magnitude to displace the degree of freedom in. Starts in the 
	            tangent plane, but can bend around edges if required. 

This function is identical to that of triangulatedMeshSpace.cpp but with different boundary condition arguments. 
*/
void fullStopMeshSpace::transportParticleAndVectors(meshPosition &pos, vector3 &displacementVector, vector<vector3> &transportVectors)
    {
    pmpFaceLocation sourceLocation = meshPositionToFaceLocation(pos);
    point3 sourcePoint = PMP::construct_point(sourceLocation,surface);
    faceIndex currentSourceFace = sourceLocation.first;

    std::vector<vertexIndex> vertexList(3);
    std::vector<point3> vertexPositions(3);

    getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
    pmpBarycentricCoordinates sourceBarycentricLocation = sourceLocation.second;
    pmpBarycentricCoordinates targetBarycentricLocation;
    bool continueShifting = true;
    vector3 currentSourceNormal = PMP::compute_face_normal(currentSourceFace,surface);
    if (abs(currentSourceNormal*displacementVector) > THRESHOLD)
        {
        displacementVector -= (currentSourceNormal*displacementVector)*currentSourceNormal;
        }
    point3 target = sourcePoint + displacementVector;
    halfedgeIndex lastUsedHalfedge(-1);
    halfedgeIndex nullHalfedge(-1);

    bool useBelowZero = true;
    
    int iter= 0;

    while(continueShifting)
        {
        iter+=1;
        //get the current barycentric coordinates of the target, clamping to the edge if they're very close 
        targetBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0],vertexPositions[1],vertexPositions[2],target);
        clampAndUpdatePosition(targetBarycentricLocation, target, currentSourceFace,surface);
        //while we don't want to eliminate negatives from the target, the source point should strictly lie within the source face, 
        //so we can use the stricter belowZeroClamp. This does not exclude the use of controls during the other parts of the loop to 
        //ensure that the source point is actually in or nearly in the face
        clampAndUpdatePosition(sourceBarycentricLocation, sourcePoint, currentSourceFace, surface,useBelowZero);
        displacementVector = vector3(sourcePoint,target);

        checkBaryNan(targetBarycentricLocation, "start of step check", iter); 
        //if the targetBarycentricLocation is in the face, we have found our destination, so check this before implementing all of the intersection locating and vector rotating logic
        bool intersectionWithEdge = false;
        for (int cc = 0; cc <3; ++cc)
            if(targetBarycentricLocation[cc] < 0)
                intersectionWithEdge = true; 
        if(!intersectionWithEdge)
            {
            continueShifting = false;
            continue;
            };
        getVertexIndicesFromFace(surface,currentSourceFace, vertexList);
        //the target barycentric location is outside the current face...find the intersection
        pmpBarycentricCoordinates intersectionPoint;
        std::vector<int> uninvolvedVertex;
        std::vector<vertexIndex> involvedVertex;

        findTriangleEdgeIntersectionInformation(sourceBarycentricLocation,targetBarycentricLocation,intersectionPoint, vertexList,lastUsedHalfedge,surface, involvedVertex,uninvolvedVertex);
        /* after above, the following are written to:
         * involvedVertex (vertex index)
         * uninvolvedVertex (integer)
         * intersection point (barycentric coords)
         */
	belowZeroClamp(intersectionPoint);                 
        point3 edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
        //walk to the location of the intersection by moving the source point. 
        vector3 toIntersection =  vector3(sourcePoint, edgeIntersectionPoint); 
	sourceBarycentricLocation = intersectionPoint;
        sourcePoint = edgeIntersectionPoint;

        if (uninvolvedVertex.size() == 2)
            {
            if (surface.is_border(vertexIndex(involvedVertex[0])))
	        {
		fullStopBoundaryVertex(sourceBarycentricLocation, targetBarycentricLocation, currentSourceNormal, currentSourceFace, displacementVector, vertexPositions, transportVectors, lastUsedHalfedge, involvedVertex[0], continueShifting);
                }
            else
                {
                updateForVertexIntersection(sourceBarycentricLocation, sourcePoint, currentSourceFace, target, displacementVector, currentSourceNormal, vertexPositions, transportVectors, lastUsedHalfedge, involvedVertex[0], toIntersection);
                continue;
                }
            }

        else if(uninvolvedVertex.size() != 1)
            {
            printSourceTargetDisplacementInfo(sourcePoint, target,sourceBarycentricLocation, targetBarycentricLocation, displacementVector);
            ERRORERROR("a barycentric coordinate of the target is negative, but neither 1 nor 2 intersections were found. Apparently some debugging is needed!");
            }

        else
            {
            halfedgeIndex intersectedEdge = surface.halfedge(involvedVertex[0], involvedVertex[1]);

            if (surface.is_border(edgeIndex(intersectedEdge)))
	        {    
		fullStopBoundaryEdge(targetBarycentricLocation, sourceBarycentricLocation, currentSourceFace, involvedVertex[0], involvedVertex[1], vertexPositions[uninvolvedVertex[0]], transportVectors, continueShifting);
                }
            else
                {
                updateForEdgeIntersection(sourceBarycentricLocation, sourcePoint, intersectionPoint, currentSourceNormal, currentSourceFace, target, transportVectors, lastUsedHalfedge, displacementVector, vertexPositions, intersectedEdge);
                continue; // probably won't be necessary in final vsn
                }

            }	
        }
    //when no more intersections, move to target (in current source face) and make sure it's legal 
    pos.faceIndex = currentSourceFace;
    pos.x = point3(targetBarycentricLocation[0],targetBarycentricLocation[1],targetBarycentricLocation[2]);   
    checkBaryNan(targetBarycentricLocation, "end of shift check" );
    }






	





