/*
Core displacement routine. Takes a degree of freedom and a direction + distance to move in, and moves 
it around the mesh, always remaining in the tangent plane of a mesh face. This particular version also transports the force on the particle along with it. 

As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex). 

pos: a meshPosition object (three doubles representing barycentric coordinates and a face index) 
     representing the position of a degree of freedom. 
displacementVector: The direction and magnitude to displace the degree of freedom in. Starts in the 
	            tangent plane, but can bend around edges if required. 
*/
void triangulatedMeshSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector, vector3 &force)
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
	clampAndUpdatePosition(targetBarycentricLocation, target, currentSourceFace);  
	//while we don't want to eliminate negatives from the target, the source point should strictly lie within the source face, 
	//so we can use the stricter belowZeroClamp. This does not exclude the use of controls during the other parts of the loop to 
	//ensure that the source point is actually in or nearly in the face
        clampAndUpdatePosition(sourceBarycentricLocation, sourcePoint, currentSourceFace, useBelowZero);
	displacementVector = vector3(sourcePoint,target);

        checkBaryNan(targetBarycentricLocation); 
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
        
	bool ifIntersection = findTriangleEdgeIntersectionInformation(sourceBarycentricLocation,targetBarycentricLocation,intersectionPoint, vertexList,lastUsedHalfedge,surface, involvedVertex,uninvolvedVertex);
        /* after above, the following are written to:
         * involvedVertex (vertex index)
         * uninvolvedVertex (integer)
         * intersection point (barycentric coords)
         */
	vector3 toIntersection;
	//below are just placeholder values which are only defined if we go through a vertex
        faceIndex provisionalTargetFace = currentSourceFace;
        vector3 newHeading = vector3(0,0,0);
        point3 edgeIntersectionPoint;
        
	if(uninvolvedVertex.size() == 2)
            {
	    cout << "Vertex intersection! " << endl;
            if (surface.is_border(vertexIndex(involvedVertex[0])))
                {
	        cout << "Boundary vertex." << endl;
                //ensure intersection point is in current face -- over the boundary does not exist.
                belowZeroClamp(intersectionPoint);                 
		edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
		//walk to the location of the intersection by moving the source point. 
		sourceBarycentricLocation = intersectionPoint;
		sourcePoint = edgeIntersectionPoint;
		    
		std::vector<faceIndex> faceIndices;
		std::vector<std::pair<vertexIndex,faceIndex>> neighborsAndFaces; 
		faceIndices.reserve(8);
                neighborsAndFaces.reserve(8); 
                
	        faceIndex cand;	
		faceCirculator fbegin(surface.halfedge(vertexIndex(involvedVertex[0])), surface), fdone(fbegin);
		do 
		    {
		    cand = *fbegin++; 
		    if (cand == currentSourceFace) continue;
		    else faceIndices.push_back(cand); 
		    } while (fbegin != fdone);       
	       
		vertexCirculator vbegin(surface.halfedge(vertexIndex(involvedVertex[0])),surface), done(vbegin);
                do 
                    {
                    for (faceIndex f: faceIndices) 
	                { 
		        if (surface.face(surface.halfedge(*vbegin++)) == f) neighborsAndFaces.push_back(make_pair(*vbegin++, f)); 
		        }	
		    } while(vbegin != done);	
		
		//below serves the rule that the next heading will be along the edge that maximally 
		//overlaps with the present heading. 
		double maxOverlap = -1; //dummy value, should always be overwritten 
		double currentOverlap = 0; 
		vector3 boundaryHeading = vector3(0,0,0); 
		faceIndex nextFace = currentSourceFace;
		
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
			}
                    }

		if (nextFace == currentSourceFace) ERRORERROR("New face not found for boundary vertex crossing.");
                    
                if (useTangentialBCs) 
		    {
		    projectOn(displacementVector,boundaryHeading);
		    target = edgeIntersectionPoint+displacementVector;	
		    }
		else 
		    {	
                    continueShifting = false;
		    targetBarycentricLocation = sourceBarycentricLocation;
		    }

		if (normalize(force)*normalize(displacementVector) > 0) projectOn(force,boundaryHeading);

		currentSourceFace = nextFace;
	        getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
                sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0], vertexPositions[1], vertexPositions[2], sourcePoint);	
	        currentSourceNormal = PMP::compute_face_normal(currentSourceFace,surface);
		
		lastUsedHalfedge = nullHalfedge;	
		checkBaryNan(sourceBarycentricLocation); //purely for debugging
                continue; //this statement ensures that this through boundary vertex branch is isolated from the below
		}
	    
	    toIntersection = vector3(sourcePoint, PMP::construct_point(std::make_pair(currentSourceFace,intersectionPoint), surface)); //always assumed to be from source to intersected vertex
            std::pair<faceIndex, vector3> targetHeading = throughVertex(involvedVertex[0], toIntersection, currentSourceFace);
            provisionalTargetFace = targetHeading.first;
            newHeading = targetHeading.second; //this heading comes out normalized
	    }

	else if(uninvolvedVertex.size() != 1)
            {
	    //below is extensive debug output in case there is a major error in shift
	    cout << endl;
	    cout << "did the program say there was  an intersection?" << endl;
	    cout << ifIntersection << endl;
	    cout << "uninvolved vertex size: " << uninvolvedVertex.size() << endl;
	    cout << "source face: " << currentSourceFace << endl;
	    cout << "source vertices: " << vertexList[0] << ", " << vertexList[1] << ", " << vertexList[2] << endl;
            printPoint(surface.point(vertexList[0]),false); 
	    printPoint(surface.point(vertexList[1]),false); 
	    printPoint(surface.point(vertexList[2]),false); 
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
	    ERRORERROR("a barycentric coordinate of the target is negative, but neither 1 nor 2 intersections were found. Apparently some debugging is needed!");
            }
	/*
        The relevant edge/vertex intersection is now identified.
        intersectionPoint contains the barycentric coordinates of the intersection point
        on a current face. Identify the next face from the relevant halfEdge, update current
        move to be from the edge intersection to the original target, and rotate it around
        the edge according to the angle of the face normals.
        */	
	        
	//update the source to be the intersection point
        vector3 targetNormal;
        halfedgeIndex intersectedEdge;
        edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
        sourceBarycentricLocation = intersectionPoint;
	sourcePoint = edgeIntersectionPoint;
        
	if (uninvolvedVertex.size() == 1)
            {
            intersectedEdge = surface.halfedge(involvedVertex[0],involvedVertex[1]);  
	    if (surface.is_border(edgeIndex(intersectedEdge)))
                {
		//intersection point *must* lie within the source face or risk a NaN at the boundary 
	        belowZeroClamp(intersectionPoint); 
	        edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
	        sourceBarycentricLocation = intersectionPoint; //we have now 'walked' to this intersection point, and will update our heading
                sourcePoint = edgeIntersectionPoint;
	        checkBaryNan(sourceBarycentricLocation);

                point3 ev1 = surface.point(involvedVertex[0]);
                point3 ev2 = surface.point(involvedVertex[1]);
		if (useTangentialBCs) 
		    { 
                    //rotate heading to be tangential to the edge in the direction most aligned with its current heading. 
		    vector3 edgeVectorForward(ev1,ev2); 
		    vector3 edgeVectorBackward(ev2,ev1); 
		    edgeVectorForward = normalize(edgeVectorForward);
	            edgeVectorBackward = normalize(edgeVectorBackward); 
		    double forwardDot = normalize(displacementVector)*edgeVectorForward;
		    double backwardDot = normalize(displacementVector)*edgeVectorBackward; 
		    double displacementLength = vectorMagnitude(displacementVector);
		    if ((forwardDot <= 0) && (backwardDot <= 0))
		         {
                         continueShifting = false;
			 } 
		    if (forwardDot > backwardDot) 
		        {
			displacementVector = forwardDot*displacementLength*edgeVectorForward; 
		        if (normalize(force)*normalize(displacementVector) > 0) projectOn(force,edgeVectorForward);
			}	
		    else 
		        {
			displacementVector = backwardDot*displacementLength*edgeVectorBackward; 
			if (normalize(force)*normalize(displacementVector) > 0) projectOn(force,edgeVectorBackward);
			}
		    //bookkeeping so as not to break later intersections
		    target = edgeIntersectionPoint+displacementVector; 
		    }
	
                //if not tangential BCs, particle stops completely if it hits a boundary
		else
		    { 
                    continueShifting = false;
		    targetBarycentricLocation = sourceBarycentricLocation; //actual displacment removal
		    vector3 orthogonalToEdge = normalize(CGAL::cross_product(currentSourceNormal, vector3(ev1,ev2)));
                    if (normalize(force)*normalize(displacementVector) > 0) projectOff(force, orthogonalToEdge);
		    }
		//source face stays the same because we can't cross boundary edges
                lastUsedHalfedge = nullHalfedge; //allows future intersections checks to include this boundary edge
		if (iter == maximumShiftEdgeCrossings) ERRORERROR("At boundary, shift proceeded for too long!"); 
		continue;
                }

	    provisionalTargetFace = surface.face(intersectedEdge);
            if (provisionalTargetFace == currentSourceFace)
                provisionalTargetFace = surface.face(surface.opposite(intersectedEdge));

            targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
            double normalDotProduct = currentSourceNormal*targetNormal;
            double angle = acos(normalDotProduct);
            if (normalDotProduct >= 1) angle = 0; //clamp for robustness against precision errors for very similar normals
	    vector3 axisVector = CGAL::cross_product(currentSourceNormal,targetNormal);
            axisVector /= vectorMagnitude(axisVector); //this could use normalize if we wanted to
            std::vector<point3> axis = {sourcePoint, sourcePoint+axisVector};
            target = rotateAboutAxis(target, axis, angle);
            displacementVector = vector3(sourcePoint, target); //move gets updated here
	    point3 forceTarget = sourcePoint + force;
	    forceTarget = rotateAboutAxis(forceTarget, axis, angle);
	    force = vector3(sourcePoint, forceTarget); 
	    //if we've gone through an edge, it can be excluded from future intersection checks
	    lastUsedHalfedge = intersectedEdge;        
       	    }

        //two conditions are exclusive
        if (uninvolvedVertex.size() == 2)
            {
            //through a vertex, we use the straightness criterion of poltier et al 2006
            double displacementLength = vectorMagnitude(displacementVector);
            displacementVector = newHeading*displacementLength;
            targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
	    target = edgeIntersectionPoint+displacementVector; 
	    //finally, set intersected edge to either of the edges containing the involved vertex so it doesn't fail to update
            intersectedEdge = surface.halfedge(involvedVertex[0],vertexList[uninvolvedVertex[0]]);
          
	    vector3 axisVector = normalize(CGAL::cross_product(currentSourceNormal,targetNormal)); 
            vector<point3> axis = {sourcePoint, sourcePoint+axisVector}; 
	    double normalDotProduct = currentSourceNormal*targetNormal;
	    double angle = acos(normalDotProduct);
	    if (normalDotProduct >= 1) angle = 0; 
	    point3 forceTarget = sourcePoint + force;
	    forceTarget = rotateAboutAxis(forceTarget, axis, angle);
	    force = vector3(sourcePoint, forceTarget);	   
	    lastUsedHalfedge = nullHalfedge; //intersecting a vertex means we need to check all edges next time
            }

        //update source face index info and new bary coords in  the new face.
        currentSourceFace = provisionalTargetFace;
        getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
        //source face has changed, so we need to update the barycentric coords of source pt
	sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0], vertexPositions[1], vertexPositions[2], sourcePoint);  	

	currentSourceNormal = targetNormal;
         
	
        if(iter == maximumShiftEdgeCrossings)
            ERRORERROR("Error: shifted across too many faces. An inadvertently large displacement was almost certainly used.");
        };

    pos.faceIndex = currentSourceFace;
    pos.x = point3(targetBarycentricLocation[0],targetBarycentricLocation[1],targetBarycentricLocation[2]);   
    checkBaryNan(targetBarycentricLocation);
    }

//Code copy-pastes a lot of the displace particle routine above...eventually refactor more nicely!
/* Big differences: 
 * -velocity bookkeeping when dealing with meshes that have boundaries is nontrivially different to ensure
 *  that velocities are updated properly
 * -velocity is rotated along with displacement as shift proceeds
 */
void triangulatedMeshSpace::transportParticleAndVelocity(meshPosition &pos, vector3 &velocityVector, vector3 &displacementVector, vector3 &force)
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
    if (abs(currentSourceNormal*velocityVector) > THRESHOLD)
        {
        velocityVector -= (currentSourceNormal*velocityVector)*currentSourceNormal;
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
	clampAndUpdatePosition(targetBarycentricLocation, target, currentSourceFace);  
	//while we don't want to eliminate negatives from the target, the source point should strictly lie within the source face, 
	//so we can use the stricter belowZeroClamp. This does not exclude the use of controls during the other parts of the loop to 
	//ensure that the source point is actually quite close to being within the face
        clampAndUpdatePosition(sourceBarycentricLocation, sourcePoint, currentSourceFace, useBelowZero);
	displacementVector = vector3(sourcePoint,target);

        checkBaryNan(targetBarycentricLocation);
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
        
	vector3 toIntersection;        
	//below are just placeholder values which are only defined if we go through a vertex
        faceIndex provisionalTargetFace = currentSourceFace;
        vector3 newHeading = vector3(0,0,0);
        point3 edgeIntersectionPoint;

	if(uninvolvedVertex.size()== 2)
            {
	    cout << "Vertex intersection!" << endl;
            if (surface.is_border(vertexIndex(involvedVertex[0])))
                {
	        cout << "Boundary vertex." << endl;
		//ensure interesction point is in current face -- over the boundary does not exist
		belowZeroClamp(intersectionPoint); 
                edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint); 
                //walk to the location of the intersection by moving the source point
		sourceBarycentricLocation = intersectionPoint;  
                sourcePoint = edgeIntersectionPoint; 

		vector<faceIndex> faceIndices;
	       	vector<pair<vertexIndex,faceIndex>> neighborsAndFaces;
                faceIndices.reserve(8); 
		neighborsAndFaces.reserve(8);

		faceIndex cand;
		faceCirculator fbegin(surface.halfedge(vertexIndex(involvedVertex[0])), surface), fdone(fbegin);
	        do 
		    {
                    cand = *fbegin++;
		    if (cand==currentSourceFace) continue;
		    else faceIndices.push_back(cand);
		    } while (fbegin != fdone);

	        vertexCirculator vbegin(surface.halfedge(vertexIndex(involvedVertex[0])),surface), vdone(vbegin);
                do 
                    {
                    for (faceIndex f: faceIndices) 
		        {
                        if (surface.face(surface.halfedge(*vbegin++)) == f) neighborsAndFaces.push_back(make_pair(*vbegin++,f));
		        }
                    } while(vbegin != vdone);
            
	    	//below serves the rule that the next heading will be along the edge that maximally 
		//overlaps with the present heading. 
		double maxOverlap = -1; //dummy value, should always be overwritten  
		double currentOverlap = 0; 
		vector3 boundaryHeading = vector3(0,0,0);
		faceIndex nextFace = currentSourceFace;

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
		       }
	           }
                
		if (nextFace == currentSourceFace) ERRORERROR("New face not found for boundary vertex crossing.");
                
		if (useTangentialBCs) 
		    {
                    projectOn(displacementVector,boundaryHeading);
		    target = edgeIntersectionPoint+displacementVector;
	            }
		else 
		    {
                    continueShifting = false;
		    targetBarycentricLocation = sourceBarycentricLocation;
		    }
                vector3 orthogonalToEdge = normalize(CGAL::cross_product(currentSourceNormal, boundaryHeading));
	        if (normalize(velocityVector)*normalize(displacementVector) > 0) projectOff(velocityVector,orthogonalToEdge); 
	        if (normalize(force)*normalize(displacementVector) > 0) projectOff(force, orthogonalToEdge); 

		currentSourceFace = nextFace; 
		getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
		sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0], vertexPositions[1], vertexPositions[2], sourcePoint);
		currentSourceNormal = PMP::compute_face_normal(currentSourceFace,surface);

		lastUsedHalfedge=nullHalfedge;
		checkBaryNan(sourceBarycentricLocation); //purely for debugging
		continue; //this statement ensures that this through boundary vertex branch is isolated from the below
		};
	    
	    //if we're not at a boundary, do standard throughVertex routine -- key function is throughVertex
            toIntersection = vector3(sourcePoint, PMP::construct_point(std::make_pair(currentSourceFace,intersectionPoint), surface)); 
            std::pair<faceIndex, vector3> targetHeading = throughVertex(involvedVertex[0], toIntersection, currentSourceFace);
            provisionalTargetFace = targetHeading.first;
            newHeading = targetHeading.second; //this heading comes out normalized
            }

	else if(uninvolvedVertex.size() != 1)
	    {
            cout << endl;
	    cout << "(within velocity transport routine)" << endl;
	    cout << "source face: " << currentSourceFace << endl;
	    cout << "source vertices: " << vertexList[0] << ", " << vertexList[1] << ", " << vertexList[2] << endl;
	    cout << surface.point(vertexList[0]) << endl; 
	    cout << surface.point(vertexList[1]) << endl; 
	    cout << surface.point(vertexList[2]) << endl; 
	    cout << "provisional target face(?): " << provisionalTargetFace << endl;
	    cout << "source r3: " << sourcePoint << endl;
            cout << "source bary: " << sourceBarycentricLocation[0] << ", " << sourceBarycentricLocation[1] << ", " << sourceBarycentricLocation[2] << endl;
            cout << "target: " << target << endl;
            cout << "target barycentric: " << targetBarycentricLocation[0] << ", " << targetBarycentricLocation[1] << ", " << targetBarycentricLocation[2] << endl;
	    cout << "displacement vector: " << displacementVector[0] << ", " << displacementVector[1] << ", " << displacementVector[2] << endl; 
	    cout << "velocity vector " << velocityVector[0] << ", " << velocityVector[1] << ", " << velocityVector[2] << endl;
	    ERRORERROR("a barycentric coordinate of the target is negative, but neither 1 nor 2 intersections were found. Apparently some debugging is needed!");
            }
        /*
        The relevant edge/vertex intersection is now identified. 
	intersectionPoint contains the barycentric coordinates of the intersection point
        on a current face. Identify the next face from the relevant halfEdge, update current
        move to be from the edge intersection to the original target, and rotate it around
        the edge according to the angle of the face normals.
        */
        	
        //update the source to be the intersection point
        vector3 targetNormal;
        halfedgeIndex intersectedEdge;
        edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
        sourceBarycentricLocation = intersectionPoint;
        sourcePoint = edgeIntersectionPoint;

	if (uninvolvedVertex.size() == 1)
            {
            intersectedEdge = surface.halfedge(involvedVertex[0],involvedVertex[1]);
            if (surface.is_border(edgeIndex(intersectedEdge)))
                {
	        //intersection point *must* lie within the source face or risk a NaN at the boundary 
	        belowZeroClamp(intersectionPoint);
	        edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
	        sourceBarycentricLocation = intersectionPoint; //we have now 'walked' to this intersection point, and will update our heading
                sourcePoint = edgeIntersectionPoint;
	        checkBaryNan(sourceBarycentricLocation);
		
	        point3 ev1 = surface.point(involvedVertex[0]);
	        point3 ev2 = surface.point(involvedVertex[1]);
		if (useTangentialBCs) 
		    { 
                    //rotate heading to be tangential to the edge in the direction most aligned with its current heading.
		    vector3 edgeVectorForward(ev1,ev2); 
		    vector3 edgeVectorBackward(ev2,ev1);
		    edgeVectorForward = normalize(edgeVectorForward);
	            edgeVectorBackward = normalize(edgeVectorBackward); 
		    double forwardDot = normalize(displacementVector)*edgeVectorForward;
		    double backwardDot = normalize(displacementVector)*edgeVectorBackward; 
		    double displacementLength = vectorMagnitude(displacementVector);
		    //if we were perpendicular to the boundary, just stop shifting
		    if ((forwardDot <= 0) && (backwardDot <= 0))
		         {
                         continueShifting = false;
			 } 
		    else if (forwardDot > backwardDot) 
		        {
			displacementVector = forwardDot*displacementLength*edgeVectorForward; 
			}	
		    else 
		        {
			displacementVector = backwardDot*displacementLength*edgeVectorBackward; 
			}
		    //bookkeeping so that we can do normal intersections later

	            //remove orthogonal components from velocity
                    vector3 orthogonalToEdge = normalize(CGAL::cross_product(currentSourceNormal, edgeVectorForward)); 
	            if (normalize(velocityVector)*normalize(displacementVector) > 0) projectOff(velocityVector,orthogonalToEdge); 
	            if (normalize(force)*normalize(displacementVector) > 0) projectOff(force, orthogonalToEdge); 
		    
		    target = edgeIntersectionPoint + displacementVector;
		    }
	
                //if not tangential BCs, particle stops completely if it hits a boundary
		else
		    {
		    //only bookkeeping necessary for non-tangential BCs is velocity -- lose all orthogonal components so we never point into non-mesh space
		    vector3 orthogonalToEdge = normalize(CGAL::cross_product(currentSourceNormal, vector3(ev1,ev2))); 
	            if (normalize(velocityVector)*normalize(displacementVector) > 0) projectOff(velocityVector,orthogonalToEdge); 
	            if (normalize(force)*normalize(displacementVector) > 0) projectOff(force, orthogonalToEdge); 
                    continueShifting = false;
		    targetBarycentricLocation = sourceBarycentricLocation;
		    }
		
		lastUsedHalfedge = nullHalfedge; 

	        if(iter ==maximumShiftEdgeCrossings ) ERRORERROR("At boundary, shift has proceeded for too long!");
		continue;
		}

	    provisionalTargetFace = surface.face(intersectedEdge);
            if (provisionalTargetFace == currentSourceFace)
                provisionalTargetFace= surface.face(surface.opposite(intersectedEdge));

            targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
            
	    double normalDotProduct = currentSourceNormal*targetNormal;
            double angle = acos(normalDotProduct);
	    if (normalDotProduct >= 1) angle = 0;

            vector3 axisVector = CGAL::cross_product(currentSourceNormal,targetNormal);
            axisVector /= vectorMagnitude(axisVector); 
            std::vector<point3> axis = {sourcePoint, sourcePoint+axisVector};
            
	    target = rotateAboutAxis(target, axis, angle);
	    displacementVector = vector3(sourcePoint, target); //move gets updated here
            
	    lastUsedHalfedge = intersectedEdge;
            //rotate the velocity vector in the same way
            point3 velocityVectorTarget = sourcePoint + velocityVector;
            velocityVectorTarget=rotateAboutAxis(velocityVectorTarget,axis,angle);
            velocityVector =  vector3(sourcePoint,velocityVectorTarget);
	    point3 forceTarget = sourcePoint + force;
	    forceTarget = rotateAboutAxis(forceTarget, axis, angle);
	    force = vector3(sourcePoint, forceTarget);
            //finally, if we've actually crossed an edge, we can eliminate it from future intersection investigation
	    lastUsedHalfedge = intersectedEdge;
            }

	//this and the condition above are exclusive
        if (uninvolvedVertex.size() == 2)
            {
            //through a vertex, we use the straightness criterion of poltier et al 2006
            double displacementLength = vectorMagnitude(displacementVector);
            displacementVector = newHeading*displacementLength;
            targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
	    target = edgeIntersectionPoint+displacementVector; 
	    //for above, we can't use a rotation with respect to the normals in question because that would violate straightness
	    //*specifically* for the path through the faces itself 
	    //for velocity vector though, it should be fine to rotate, as rotation is how we transport vector quantity into any new face
	    vector3 axisVector = CGAL::cross_product(currentSourceNormal,targetNormal);
	    axisVector /= vectorMagnitude(axisVector); // normalize works too
            std::vector<point3> axis = {sourcePoint, sourcePoint+axisVector};
	    double normalDotProduct = currentSourceNormal*targetNormal;
            double angle = acos(normalDotProduct);
	    if (normalDotProduct >= 1) angle = 0;
	    point3 velocityVectorTarget = sourcePoint + velocityVector;
            velocityVectorTarget=rotateAboutAxis(velocityVectorTarget,axis,angle);
            velocityVector =  vector3(sourcePoint,velocityVectorTarget);
	    point3 forceTarget = sourcePoint + force;
	    forceTarget = rotateAboutAxis(forceTarget, axis, angle);
	    force = vector3(sourcePoint, forceTarget);
	    //we didn't intersect a halfedge, so we can't exclude one from future movements
	    lastUsedHalfedge = nullHalfedge; 
            }

        //update source face index info and new bary coords in  the new face.
        currentSourceFace = provisionalTargetFace; //we now know we've moved into this face, and case use it as the source
        
	getVertexPositionsFromFace(surface, currentSourceFace, vertexPositions);
        //source face has changed, so we need to update the barycentric coords of source pt
	sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0], vertexPositions[1], vertexPositions[2], sourcePoint);  	

        currentSourceNormal = targetNormal;
        
	if(iter ==maximumShiftEdgeCrossings ) ERRORERROR("Error: shifted across an extremely  large number of faces... is this an error, or an accidental call with an extremely large displacement vector?");
        };

    pos.faceIndex = currentSourceFace;
    pos.x = point3(targetBarycentricLocation[0],targetBarycentricLocation[1],targetBarycentricLocation[2]);
    };

