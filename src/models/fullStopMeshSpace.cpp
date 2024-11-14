void fullStopMeshSpace::printSourceTargetDiplacementInfo(point3 sourcePoint, point3 target, pmpBarycentricCoordinates sourceBarycentricLocation, pmpBarycentricCoordinates targetBarycentricLocation, vector3 displacementVector)
            {
            //below is extensive debug output in case there is a major error in shift
            cout << endl;
            cout << "source r3: "; 
            printPoint(sourcePoint); 
            cout << endl;

            cout << "source bary: ";
            printBary(sourceBarycentricLocation, true); 
            cout << endl;

            cout << "target: "; 
            printPoint(target, true); 
            cout << endl;

            cout << "target barycentric: ";
            printBary(targetBarycentricLocation,true);
            cout << endl;

            cout << "displacement vector: " << displacementVector[0] << ", " << displacementVector[1] << ", " << displacementVector[2] << endl; 
            }


//below also updates source face -- perhaps I can figure out a better name? 
pair<vector3, vertexIndex> fullStopMeshSpace::getBoundaryVertexHeading(vertexIndex v, faceIndex &sourceFace, pmpBarycentricCoordinates &sourceBarycentricLocation, vector<point3> vertexPositions, vector3 sourceNormal) 
    {
    std::vector<faceIndex> faceIndices;
    std::vector<std::pair<vertexIndex,faceIndex>> neighborsAndFaces; 
    faceIndices.reserve(8);
    neighborsAndFaces.reserve(8); 

    faceIndex cand;	
    faceCirculator fbegin(surface.halfedge(vertexIndex(involvedVertex[0])), surface), fdone(fbegin);
    do 
        {
        cand = *fbegin++; 
        if (cand == currentSourceFace) 
            continue;
        else 
            faceIndices.push_back(cand); 
        } while (fbegin != fdone);       

    vertexCirculator vbegin(surface.halfedge(vertexIndex(involvedVertex[0])),surface), done(vbegin);
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
    faceIndex nextFace = currentSourceFace;
    vertexIndex boundaryVertex = involvedVertex[0]; 

    for (auto vi: neighborsAndFaces) 
        {
        vector3 outwardVector(surface.point(involvedVertex[0]),surface.point(vi.first));
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

    if (nextFace == currentSourceFace)
        ERRORERROR("New face not found for boundary vertex crossing.");
    if (boundaryVertex == involvedVertex[0]) 
        ERRORERROR("vertex not updated in throughvertex boundary crossing");
   
    // last step is doing small updates to source face, which happen by address here 
    // as it is the most convenient place to do so -- technically none of this is necessary
    // for full stop  
    sourceFace = nextFace;
    getVertexPositionsFromFace(surface, sourceFace, vertexPositions);
    sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0], vertexPositions[1], vertexPositions[2], sourcePoint);
    sourceNormal = PMP::compute_face_normal(sourceFace, surface);    
    return make_pair(boundaryHeading, boundaryVertex);
    }

void fullStopMeshSpace::projectVectorsForBoundaryVertexIntersection(vector3 heading, vertexIndex interesctedV, vertexIndex boundaryV, vector3 fNormal, faceIndex face, vector<vector3> &transportVectors) 
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
        vector3 inwardVector(surface.point(involvedVertex[0]), insideVertex); 
        projectVectorsIfOverBoundary(transportVectors, orthogonalToEdge, inwardVector);	
    }


void fullStopMeshSpace::updateForBoundaryVertexIntersection(pmpBarycentricCoordinates &sourceBarycentricLocation, point3 &sourcePoint, vector3 &sourceNormal, pmpBarycentricCoordinates &intersectionPoint, faceIndex &currentSourceFace, point3 &target, vector<point3> &vertexPositions, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vertexIndex intersectedV);
    {
    //first get heading & next vertex along the boundary in that heading; source face is modified
    pair<vector3, vertexIndex> headingAndVertex = getBoundaryVertexHeading(intersectedV, currentSourceFace, sourceBarycentricLocation, vertexPositions, sourceNormal);
    //project any vectors we need to project according to the calculated boundary heading
    if (transportVectors.size() > 0)
        projectVectorsForBoundaryVertexIntersection(headingAndVertex.first, intersectedV, headingAndVertex.second, sourceNormal, currentSourceFace, transportVectors);
    
    //technically, the above moves us into an adjacent face -- since the vertex is in the same place
    //either way, it doesn't change the outcome, but this version is more extensible to tangential
    //BCs at the cost of more code overhead.  
    continueShifting = false; 
    targetBarycentricLocation = sourceBarycentricLocation;
     
    halfedgeIndex nullHalfedge(-1); 
    lastUsedHalfedge = nullHalfedge;  
    }


void fullStopMeshSpace::updateForVertexIntersection(pmpBarycentricCoordinates &sourceBCs, point3 &sourcePoint, faceIndex &sourceFace, point3 &target, vector3 &displacementVector, vector3 &sourceNormal, vector<point3> vertexPositions, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vertexIndex intersectedV, vector3 toIntersection)
    {
    //through vertex call to enable updates
    pair<faceIndex, vector3> targetHeading = throughVertex(intersectedV, toIntersection, sourceFace);
    faceIndex targetFace = targetHeading.first; 
    vector3 newHeading = targetHeading.second; 
   
    //main updates
    double displacementLength = vectorMagnitude(displacementVector);
    vector3 targetNormal = PMP::compute_face_normal(targetFace,surface);
    displacementVector = newHeading*displacementLength;
    target = sourcePoint+displacementVector; 

    //we don't need any of these if there are no vectors to transport
    if (transportVectors.size() > 0)
        {  
        vector3 axisVector = normalize(CGAL::cross_product(sourceNormal,targetNormal)); 
        vector<point3> axis = {sourcePoint, sourcePoint+axisVector}; 
        double normalDotProduct = currentSourceNormal*targetNormal;
        double angle = acos(normalDotProduct);
        if (normalDotProduct >= 1) angle = 0; 

        for (int i = 0; i < transportVectors.size(); i++)
            {
            vector3 vTr = transportVectors[i];
            point3 vTrTarget = sourcePoint + vTr;
            vTrTarget = rotateAboutAxis(vTrTarget, axis, angle);
            transportVectors[i] = vector3(sourcePoint, vTrTarget);
            }
        }	   
    //update source face index info and new bary coords in  the new face.
    sourceFace = targetFace;
    getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
    //source face has changed, so we need to update the barycentric coords of source pt
    sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0], vertexPositions[1], vertexPositions[2], sourcePoint);  	

    sourceNormal = targetNormal;
    lastUsedHalfedge = nullHalfedge;    
    }




//"easy" one -- just rotate over edge and move target. 
void fullStopMeshSpace::updateForEdgeIntersection(pmpBarycentricCoordinates &sourceBarycentricLocation, point3 &sourcePoint, pmpBarycentricCoordinates &intersectionPoint, vector3 &currentSourceNormal, faceIndex &currentSourceFace, point3 &target, vector<vector3> &transportVectors, halfedgeIndex &lastUsedHalfedge, vector<point3> &vertexPositions, halfedgeIndex intersectedEdge)
    {
    sourcePoint = globalSMSP->point(currentSourceFace, intersectionPoint);

    faceIndex provisionalTargetFace = surface.face(intersectedEdge);
    if (provisionalTargetFace == currentSourceFace)
        provisionalTargetFace = surface.face(surface.opposite(intersectedEdge));

    vector3 targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
    double normalDotProduct = currentSourceNormal*targetNormal;
    double angle = acos(normalDotProduct);
    if (normalDotProduct >= 1) 
        angle = 0; //clamp for robustness against precision errors for very similar normals
    vector3 axisVector = CGAL::cross_product(currentSourceNormal,targetNormal);
    axisVector /= vectorMagnitude(axisVector); //this could use normalize if we wanted to
    std::vector<point3> axis = {sourcePoint, sourcePoint+axisVector};

    target = rotateAboutAxis(target, axis, angle);
    displacementVector = vector3(sourcePoint, target); //move gets updated here

    for (int i = 0; i < transportVectors.size(); i++) 
        {
        vector3 vTr = transportVectors[i];
        point3 vTrTarget = sourcePoint + vTr;
        vTrTarget = rotateAboutAxis(vTrTarget, axis, angle);
        transportVectors[i] = vector3(sourcePoint, vTrTarget);
        }

    //if we've gone through an edge, it can be excluded from future intersection checks
    lastUsedHalfedge = intersectedEdge;

    //update source face index info and new bary coords in  the new face.
    currentSourceFace = provisionalTargetFace;
    getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
    //source face has changed, so we need to update the barycentric coords of source pt
    sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0], vertexPositions[1], vertexPositions[2], sourcePoint);  	

    currentSourceNormal = targetNormal;
    }

void fullStopMeshSpace::fullStopAtBoundaryEdge(pmpBarycentricCoordinates &targetBCs, pmpBarycentricCoordinates &sourceBCs, bool continueShifting, vertexIndex edgeV1, vertexIndex edgeV2, point3 innerVertex, transportVectors);
    {
    if (transportVectors.size() > 0) 
        {
	//below handles just tangential BCs and/or nonzero transport vectors size -- later
	//iterations will avoid always executing it 
        point3 ev1 = surface.point(edgeV1);
        point3 ev2 = surface.point(edgeV2);
                     
        //rotate heading to be tangential to edge in the direction most aligned with current heading 
        vector3 edgeVectorForward(ev1,ev2); 
        vector3 edgeVectorBackward(ev2,ev1); 
        edgeVectorForward = normalize(edgeVectorForward);
        edgeVectorBackward = normalize(edgeVectorBackward); 
        double forwardDot = normalize(displacementVector)*edgeVectorForward;
        double backwardDot = normalize(displacementVector)*edgeVectorBackward; 
        double displacementLength = vectorMagnitude(displacementVector);

        //check to see if any of our transported vectors were pointing over a boundary along with the displacement & project if so 
        vector3 orthogonalToEdge = normalize(CGAL::cross_product(currentSourceNormal, vector3(ev1,ev2)));
        //if orthogonal is pointing out, then we want to project only if the overlap with the orthogonal vector is positive; 
        //if orthogonal is pointing in, then we want to project only if the overlap with the orthogonal vector is negative. 
        vector3 inwardVector(edgeIntersectionPoint, innerVertex); 
        projectVectorsIfOverBoundary(transportVectors, orthogonalToEdge, inwardVector);
	}	
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
        edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
        //walk to the location of the intersection by moving the source point. 
        vector3 toIntersection =  vector3(sourcePoint, edgeIntersectionPoint); 
	sourceBarycentricLocation = intersectionPoint;
        sourcePoint = edgeIntersectionPoint;

	if (uninvolvedVertex.size() == 2)
	    {
	    if (surface.is_border(vertexIndex(involvedVertex[0]) i
	        {
	        updateForBoundaryVertexIntersection(sourceBarycentricLocation, sourcePoint, intersectionPoint, currentSourceFace, target, transportVectors, lastUsedHalfedge);
		continue;
	        }
	    else 
	        {
	        updateForVertexIntersection(sourceBarycentricLocation, sourcePoint, currentSourceFace, target, intersectionPoint, lastUsedHalfedge);
		continue;
		}
	    }
        
	else if(uninvolvedVertex.size() != 1)
            {
	    printSourceTargetDisplacementInfo(sourcePoint, sourceBarycentricLocation, target, targetBarycentricLocation, displacementVector); 
            ERRORERROR("a barycentric coordinate of the target is negative, but neither 1 nor 2 intersections were found. Apparently some debugging is needed!");
            }

	else 
	    {
	    intersectedEdge = surface.halfedge(involvedVertex[0], involvedVertex[1]); 
            if (surface.is_border(edgeIndex(intersectedEdge)))
                {
                fullStopAtBoundaryEdge(targetBarycentricLocation, sourceBarycentricLocation, continueShifting, involvedVertex[0], involvedVertex[1], vertexPositions[uninvolvedVertex[0]], transportVectors);
                continue;        
		}
	    else 
	        {
                updateForEdgeIntersection(sourceBarycentricLocation, sourcePoint, intersectionPoint, transportVectors, lastUsedHalfedge, currentSourceFace, intersectedEdge);
		continue; // probably won't be necessary in final vsn
                }	

	    }
        }
    //when no more intersections, move to target (in current source face) and make sure it's legal 
    pos.faceIndex = currentSourceFace;
    pos.x = point3(targetBarycentricLocation[0],targetBarycentricLocation[1],targetBarycentricLocation[2]);   
    checkBaryNan(targetBarycentricLocation, "end of shift check" );
    }






	





