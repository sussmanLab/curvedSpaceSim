#include "triangulatedMeshSpace.h"

/*! \file triangulatedMeshSpace.cpp */
#include <stdexcept>

void triangulatedMeshSpace::updateMeshSpanAndTree()
    {
    //set domain in which surface lives
    minVertexPosition.x = 0;minVertexPosition.y = 0;minVertexPosition.z = 0;
    maxVertexPosition.x = 0;maxVertexPosition.y = 0;maxVertexPosition.z = 0;
    for (vertexIndex v : surface.vertices())
        {
        point3 p = surface.point(v);                          
	if(p[0] < minVertexPosition.x)
            minVertexPosition.x = p[0];
        if(p[1] < minVertexPosition.y)
            minVertexPosition.y = p[1];
        if(p[2] < minVertexPosition.z)
            minVertexPosition.z = p[2];
        if(p[0] > maxVertexPosition.x)
            maxVertexPosition.x = p[0];
        if(p[1] > maxVertexPosition.y)
            maxVertexPosition.y = p[1];
        if(p[2] > maxVertexPosition.z)
            maxVertexPosition.z = p[2];
        };
    if(verbose)
        printf("mesh spans (%f,%f,%f) to (%f,%f,%f)\n", minVertexPosition.x,minVertexPosition.y,minVertexPosition.z, maxVertexPosition.x,maxVertexPosition.y,maxVertexPosition.z);

    globalSMSP = make_shared<surfaceMeshShortestPath>(surface);
    //global tree speeds up all subsequent location operations on the surface
    globalSMSP->build_aabb_tree(globalTree);
    }

void triangulatedMeshSpace::loadMeshFromFile(std::string filename, bool _verbose)
    {
    verbose = _verbose;
    positionsAreEuclidean = false;
    surface = triangleMesh(); //clear mesh so we don't accidentally load two meshes at once
    if(verbose)
        {
        printf("loading from file %s\n",filename.c_str());
        }
    if(!CGAL::IO::read_polygon_mesh(filename, surface))
        {
        if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, surface))
            {
            std::cerr << "Invalid input file." << std::endl;
            throw std::exception();
            }
        };
    if(!CGAL::is_triangle_mesh(surface))
        {
        std::cerr << "Non-triangular mesh" << std::endl;
        throw std::exception();
        };
    int nFaces = surface.number_of_faces();
    int nVertices = surface.number_of_vertices();
    if(verbose)
        {
        printf("input mesh has %i faces and %i vertices\n",nFaces,nVertices);
        };
    //set spatial domain for surface, create global AABB tree
    updateMeshSpanAndTree();
    };

void triangulatedMeshSpace::isotropicallyRemeshSurface(double targetEdgeLength)
    {
    cout << "number of vertices before remesh: " << surface.number_of_vertices() << endl;
    PMP::isotropic_remeshing(surface.faces(), targetEdgeLength, surface);
    cout << "number of vertices after  remesh: " << surface.number_of_vertices() << endl;
    updateMeshSpanAndTree();
    };

void triangulatedMeshSpace::meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result)
    {
    if(result.size()!=p1.size())
        result.resize(p1.size());
    for(int ii = 0; ii < p1.size();++ii)
        {
        smspFaceLocation sourcePoint = meshPositionToFaceLocation(p1[ii]);
        point3 b=globalSMSP->point(sourcePoint.first,sourcePoint.second);
        result[ii].x=b;
        };
    }

void triangulatedMeshSpace::meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result)
    {
    if(result.size()!=p1.size())
        result.resize(p1.size());
    for(int ii = 0; ii < p1.size();++ii)
        {
        smspFaceLocation sourcePoint = meshPositionToFaceLocation(p1[ii]);
        point3 b=globalSMSP->point(sourcePoint.first,sourcePoint.second);
        result[ii].x=b[0];
        result[ii].y=b[1];
        result[ii].z=b[2];
        };
    }

void triangulatedMeshSpace::randomPosition(meshPosition &p, noiseSource &noise)
    {
    double3 baryPoint;
    baryPoint.x=noise.getRealUniform();
    baryPoint.y=noise.getRealUniform(0,1-baryPoint.x);
    baryPoint.z=1-baryPoint.x-baryPoint.y;
    p.x = point3(baryPoint.x,baryPoint.y,baryPoint.z);
    p.faceIndex = noise.getInt(0,surface.number_of_faces()-1);
    }

void triangulatedMeshSpace::randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise)
    {
    //compute the normal to the face containing point p
    pmpFaceLocation sourceLocation = meshPositionToFaceLocation(p);
    faceIndex sourceFace = sourceLocation.first;
    vector3 sourceNormal = PMP::compute_face_normal(sourceFace,surface);

    //construct an orthogonal vector by hand
    vector3 orthogonalizer;
    if(sourceNormal[0] ==0 && sourceNormal[1] ==0 )
        orthogonalizer = vector3(sourceNormal[1]-sourceNormal[2], sourceNormal[2]-sourceNormal[0], sourceNormal[0]-sourceNormal[1]);
    else
        orthogonalizer = vector3(sourceNormal[1],-sourceNormal[0],0);
    //make this orthogonal vector a unit vector
    orthogonalizer /= sqrt(orthogonalizer.squared_length());
    //grab the second tangent vector (which should already be a unit vector)
    vector3 tangent2 =CGAL::cross_product(sourceNormal,orthogonalizer);

    //finally, choose a gaussian weight of each of the two orthogonal vectors in the face's tangent
    v = noise.getRealNormal()*orthogonalizer+noise.getRealNormal()*tangent2;
    }

void triangulatedMeshSpace::convertToEuclideanPositions(std::vector<meshPosition> &a, std::vector<meshPosition> &b)
    {
    int N = a.size();
    if(b.size()!=N)
        b.resize(N);
    for(int ii = 0; ii < N; ++ii)
        {
        smspFaceLocation sourcePoint = meshPositionToFaceLocation(a[ii]);
        b[ii].x=globalSMSP->point(sourcePoint.first,sourcePoint.second);
        b[ii].faceIndex = a[ii].faceIndex;
        };
    }


void triangulatedMeshSpace::distanceWithSubmeshing(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent,double distanceThreshold)
    {
    //create submesh:
    smspFaceLocation sourcePoint = meshPositionToFaceLocation(p1);
    std::vector<smspFaceLocation> faceTargetsForSubmesh(p2.size());
    for(int ii = 0; ii < p2.size(); ++ii)
        faceTargetsForSubmesh[ii] = meshPositionToFaceLocation(p2[ii]);
    std::unordered_map<faceIndex,int> faceMap;
    std::unordered_map<vertexIndex,int> vertexMap;
    double currentDistanceThreshold  = maximumDistance;
    if(distanceThreshold<maximumDistance)
        currentDistanceThreshold = distanceThreshold;

    triangleMesh submesh = submeshAssistant.constructSubmeshFromSourceAndTargets(surface, sourcePoint,faceTargetsForSubmesh,currentDistanceThreshold,vertexMap,faceMap);

    int nTargets = p2.size();
    distances.resize(nTargets);
    startPathTangent.resize(nTargets);
    endPathTangent.resize(nTargets);

    //Note that in the process of submeshing, the original barycentric coordinates may get permuted
    //in the submesh (i.e., the ordering of vertices around a face is not guaranteed to be preserved)
    //and the faceIndex needs to change
    convertBarycentricCoordinates(surface,submesh,faceMap,sourcePoint);
    for(int ii = 0; ii < nTargets; ++ii)
        {
        convertBarycentricCoordinates(surface,submesh,faceMap,faceTargetsForSubmesh[ii]);
        }
    //now that indexing is all re-aligned, create local surfaceMeshShortestPath
    surfaceMeshShortestPath localSMSP(submesh);
    AABB_tree localTree;
    localSMSP.build_aabb_tree(localTree);
    localSMSP.add_source_point(sourcePoint.first,sourcePoint.second);
    localSMSP.build_sequence_tree();
    for(int ii = 0; ii < nTargets; ++ii)
        {
        computePathDistanceAndTangents(&localSMSP, faceTargetsForSubmesh[ii], distances[ii],startPathTangent[ii], endPathTangent[ii]);
        //if the submesh has multiple connected components, the distance will be returned as negative
        if(distances[ii] <0)
            {
            distances[ii] = 2.0*maximumDistance;
            startPathTangent[ii] = {0,0,1};
            endPathTangent[ii] = {0,0,1};
            /*
            printf("disconnected submesh %f %f \n",maximumDistance,distanceThreshold);//for debugging more aggresive submeshing
            printf("vertices={");
            for(vertexIndex v : submesh.vertices())
                {
                    printPoint(submesh.point(v));printf(",");
                }
            printf("};\n");
            printf("target={");
            printPoint(PMP::construct_point(faceTargetsForSubmesh[ii],submesh)); printf("};\n");
            printf("source={"); printPoint(PMP::construct_point(sourcePoint,submesh)); printf("};\n");
            */
            }
        };
    };

/*
As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex)
*/
void triangulatedMeshSpace::distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent, double distanceThreshold)
    {
    if(submeshingActivated)
        {
        distanceWithSubmeshing(p1,p2,distances,startPathTangent,endPathTangent,distanceThreshold);
        return;
        }
    int nTargets = p2.size();
    distances.resize(nTargets);
    startPathTangent.resize(nTargets);
    endPathTangent.resize(nTargets);

    smspFaceLocation sourcePoint = meshPositionToFaceLocation(p1);

    globalSMSP->add_source_point(sourcePoint.first,sourcePoint.second);
    globalSMSP->build_sequence_tree();

    for(int ii = 0; ii < nTargets; ++ii)
        {
        smspFaceLocation targetPoint = meshPositionToFaceLocation(p2[ii]);
        computePathDistanceAndTangents(globalSMSP.get(), targetPoint, distances[ii],startPathTangent[ii], endPathTangent[ii]);
        };
    globalSMSP->remove_all_source_points();
    };

/*!throughVertex takes as input the index of the vertex we're going through, the vector from the current source to that intersection,
and the source face. Then, it cycles through the faces adjoining the vertex, collecting the total angle of the vertex (the total
angle subtended by the edge pairs), and determines the heading that would be halfway through the total angle. The face containing
that heading and the heading itself are returned. The half-of-total angle criterion is the same as the straightest geodesic criterion;
the straightest path has equal angles on both of its sides (e.g. in flat space, 180 degrees either way you rotate).  
*/
std::pair<faceIndex,vector3> triangulatedMeshSpace::throughVertex(vertexIndex &intersectedVertex, vector3 &toIntersection, faceIndex &sourceFace)  
    {
    point3 vi_r3 = surface.point(intersectedVertex);

    //determine what vertex should be the starting vertex for the angle measurement as
    //the next vertex from the intersected one (guaranteeing it's in the original face)
    vertexIndex startingVertex;

    //if we know what the source face is, this guarantees the first edge
    //we consider will be part of the source face.
    halfedgeIndex hf = surface.halfedge(sourceFace);
    for(halfedgeIndex hi : halfedges_around_face(hf, surface))
        {
        if(source(hi, surface) == intersectedVertex) startingVertex = target(hi,surface);
        }

    std::vector<vertexIndex> neighborIndices;
    neighborIndices.reserve(8); //this should be as large as is possible -- on average, expect 6

    vertexCirculator vbegin(surface.halfedge(intersectedVertex),surface), done(vbegin);
    do 
        {
        neighborIndices.push_back(*vbegin++);
        } while(vbegin != done);

    std::vector<point3> neighborVertices;
    std::vector<vertexIndex> neighborIndicesCopy; //going to sort this one by starting point first
    neighborVertices.reserve(neighborIndices.size());
    neighborIndicesCopy.reserve(neighborIndices.size());
    vertexIndex currentVert;
    bool started = false;
    bool atStartPoint = false;
    bool unfinished = true;
    int ii = -1;
    //we want a list that retains the original ordering, but 
    //has a new start point in a totally mutable way (e.g. without
    //being able to guess where the start point is in the original list)
    while (unfinished) 
        {
        //ii used to iterate through list
        ii = (ii+1)%neighborIndices.size(); //we need to be able to wrap if the start point is late in the list
        currentVert = neighborIndices[ii]; 
        atStartPoint = (currentVert == startingVertex); //if we're at the right starting point...
        if (started && atStartPoint) unfinished=false; //finishing counting, or
        else //start counting + continue counting 
       	    {
            if (atStartPoint) started=true; //this lets us know we can start counting, second condition will be filled until done
            if (started)
	        {
                neighborIndicesCopy.push_back(currentVert);
                neighborVertices.push_back(surface.point(currentVert));
                };
            };
        }

    std::vector<vector3> edgeVectors;
    int numNeighbors = neighborVertices.size();
    edgeVectors.reserve(numNeighbors);

    for (int i = 0; i < numNeighbors; i++) {
        point3 v = neighborVertices[i];
        edgeVectors.push_back(normalize(vector3(vi_r3, v)));
    }

    double totalAngle = 0;
    for (int i = 0; i < numNeighbors; i++)
        {
        vector3 v1 = edgeVectors[i];
        vector3 v2 = edgeVectors[(i+1)%numNeighbors];
        totalAngle += angleBetween(v1,v2);
        }
     //now -- pick the vector in a target face that corresponds to having
    //exactly half of the total angle of the vertex between it and the
    //vector from the original source to the vertex itself.
    double angleToTravel = totalAngle/2.0;
    vector3 queryHeading = normalize(-toIntersection); // negative so it matches heading of edge vectors
    double firstAngle = angleBetween(queryHeading,edgeVectors[0]);
    double angleTraveled = firstAngle;
    int counter = 0;
    // below always terminates because of final counter check
    while (angleTraveled < angleToTravel) 
        {
        vector3 e1 = edgeVectors[counter];
        vector3 e2 = edgeVectors[(counter+1)%numNeighbors]; //leaving this here for instances where one face is gigantic and the others are miniscule
        angleTraveled += angleBetween(e1, e2);
        counter ++; // counter will always finish on the e2 index
        if (counter >= numNeighbors) 
	    {
            printf("did not find an angular halfway when picking heading from a vertex, debug!");
            break;
            }
        }
    double angleDifference = angleToTravel - angleTraveled;
    vector3 rotAxisVec = normalize(CGAL::cross_product(edgeVectors[counter-1], edgeVectors[counter]));
    std::vector<point3> rotAxis;
    rotAxis.reserve(2);
    rotAxis.push_back(vi_r3);
    rotAxis.push_back(vi_r3 + rotAxisVec);
    point3 rotationTarget = vi_r3 + edgeVectors[counter];

    point3 prospPoint = rotateAboutAxis(rotationTarget, rotAxis, angleDifference);
    vector3 prospHeading = normalize(vector3(vi_r3,prospPoint));
    
    //now -- need to select target face as face subtended by e1 and e2
    std::vector<faceIndex> candidateFaces;
    candidateFaces.reserve(8); // should be six, but managing pathological cases
    faceCirculator fbegin(surface.halfedge(intersectedVertex),surface), fdone(fbegin);
    faceIndex cand;
    do 
        {
        cand = *fbegin++; //the next face to add is the face pointed to by the next
                        //iteration of the iterator
        if (cand == sourceFace) continue; //don't need to consider the face we're coming from
        else 
	    {
            candidateFaces.push_back(cand);
            }
        } while(fbegin != fdone);

    //null result is target face being the same as the source face
    faceIndex targetFace = sourceFace;
    std::set<vertexIndex> targetFaceVerts;
    targetFaceVerts.insert(intersectedVertex);
    targetFaceVerts.insert(neighborIndicesCopy[counter-1]);
    targetFaceVerts.insert(neighborIndicesCopy[counter]);

    for (faceIndex face: candidateFaces) 
        {
        //below is glorified check to see if all the vertices of the face are the 
	//same as the vertices of the two edges we need to be included
        std::vector<vertexIndex> currentFaceVs(3); 
	getVertexIndicesFromFace(surface, face, currentFaceVs);
        std::set<vertexIndex> currentVsSet;
        currentVsSet.insert(currentFaceVs[0]);
        currentVsSet.insert(currentFaceVs[1]);
        currentVsSet.insert(currentFaceVs[2]);
        bool allThere = true;
        for (vertexIndex vert: targetFaceVerts)
            {
            auto it = currentVsSet.find(vert);
            if (it == currentVsSet.end())
                {
                allThere = false;
                break;
                };
            };
        if (allThere == true)
            {
            targetFace=face;
            break;
            }
        }
    if (sourceFace == targetFace) 
        {
	printf("throughvertex didn't find a face, debug!"); 
	throw std::exception();
        }
    return std::make_pair(targetFace,prospHeading); 
    }  


/*
Core displacement routine. Takes a degree of freedom and a direction + distance to move in, and moves 
it around the mesh, always remaining in the tangent plane of a mesh face.

As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex). 

pos: a meshPosition object (three doubles representing barycentric coordinates and a face index) 
     representing the position of a degree of freedom. 
displacementVector: The direction and magnitude to displace the degree of freedom in. Starts in the 
	            tangent plane, but can bend around edges if required. 
*/
void triangulatedMeshSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector)
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
        //printf("%g,\n", abs(currentSourceNormal*displacementVector));
        displacementVector -= (currentSourceNormal*displacementVector)*currentSourceNormal;
        //printf("%g,\n", abs(currentSourceNormal*displacementVector));
        //ERRORERROR("non-tangent displacement vector on a face");
        }
    point3 target = sourcePoint + displacementVector;
    vector3 currentMove = vector3(sourcePoint, target);
    halfedgeIndex lastUsedHalfedge(-1);

int iter= 0;
    while(continueShifting)
        {
iter+=1;
        //get the current barycentric coordinates of the target
        targetBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0],vertexPositions[1],vertexPositions[2],target);
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
	vector3 toIntersection;
        //below are just placeholder values which are only defined if we go through a vertex
        faceIndex provisionalTargetFace = currentSourceFace;
        vector3 newHeading = vector3(0,0,0);

        if(uninvolvedVertex.size()== 2)
            {
	    //if we're at a border, we need to invoke boundary conditions -- for now, this is just
	    //stopping the particle cold. is_border works on both edges and vertices
            if (surface.is_border(vertexIndex(involvedVertex[0])))
                {
                targetBarycentricLocation = intersectionPoint;
                //clamp so we are not actually on the edge (or over it!), but just very close -- avoids 
                //confusion on the intersection checker's part at later timesteps. 
                clampToThreshold(targetBarycentricLocation);
                continueShifting = false;
                continue;
                }
	    toIntersection = vector3(sourcePoint, PMP::construct_point(std::make_pair(currentSourceFace,intersectionPoint), surface)); //always assumed to be from source to intersected vertex
            printf("updating shift direction with throughvertex\n");
            std::cout << "involved vertex " << involvedVertex[0] << std::endl;
            std::pair<faceIndex, vector3> targetHeading = throughVertex(involvedVertex[0], toIntersection, currentSourceFace);
            provisionalTargetFace = targetHeading.first;
            newHeading = targetHeading.second; //this heading comes out normalized
            }
	if(uninvolvedVertex.size() != 1)
            ERRORERROR("a barycentric coordinate of the target is negative, but neither 1 nor 2 intersections were found. Apparently some debugging is needed!");
        /*
        We have now identified the relevant edge or vertex intersection point.
        intersectionPoint contains the barycentric coordinates of the intersection point
        on a current face. Identify the next face from the relevant halfEdge, update current
        move to be from the edge intersection to the original target, and rotate it around
        the edge according to the angle of the face normals.
        */
        point3 edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
        sourcePoint = edgeIntersectionPoint;

        currentMove = vector3(edgeIntersectionPoint, target);

        //identify the next faceIndex from the shared intersected edge
        vector3 targetNormal;
        halfedgeIndex intersectedEdge;
        if (uninvolvedVertex.size() == 1)
            {
            intersectedEdge = surface.halfedge(involvedVertex[0],involvedVertex[1]);  
            //same as above -- particle stops completely if it hits a boundary, for now
	    if (surface.is_border(edgeIndex(intersectedEdge)))
                {
                targetBarycentricLocation = intersectionPoint;
                clampToThreshold(targetBarycentricLocation);
                continueShifting = false;
                continue;
                }

	    provisionalTargetFace = surface.face(intersectedEdge);
            if (provisionalTargetFace == currentSourceFace)
                provisionalTargetFace= surface.face(surface.opposite(intersectedEdge));

            targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
            double normalDotProduct = currentSourceNormal*targetNormal;
            double angle = acos(normalDotProduct);
            if (normalDotProduct >= 1) angle = 0; //clamp for robustness against precision errors for very similar normals
	    vector3 axisVector = CGAL::cross_product(currentSourceNormal,targetNormal);
            axisVector /= vectorMagnitude(axisVector); //this could use normalize if we wanted to
            std::vector<point3> axis = {sourcePoint, sourcePoint+axisVector};
            target = rotateAboutAxis(target, axis, angle);
            displacementVector = vector3(sourcePoint, target); //move gets updated here
            //we can only update last used halfedge when crossing an edge -- maybe ask daniel if this breaks for going through vertices
            }

        //two conditions are exclusive
        if (uninvolvedVertex.size() == 2)
            {
            //through a vertex, we use the straightness criterion of poltier et al 2006
            double displacementLength = vectorMagnitude(displacementVector);
            displacementVector = newHeading*displacementLength;
            targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
	    //finally, set intersected edge to either of the edges containing the involved vertex so it doesn't fail to update
            intersectedEdge = surface.halfedge(involvedVertex[0],vertexList[uninvolvedVertex[0]]);
            }

        //update source face index info and new bary coords in  the new face.
        currentSourceFace = provisionalTargetFace;
        getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
        sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0],vertexPositions[1],vertexPositions[2],edgeIntersectionPoint);
        currentSourceNormal = targetNormal;
        lastUsedHalfedge = intersectedEdge;        
	
        if(iter == maximumShiftEdgeCrossings)
            ERRORERROR("Error: shifted across too many faces. An inadvertently large displacement was almost certainly used.");
        };

    pos.faceIndex = currentSourceFace;
    pos.x = point3(targetBarycentricLocation[0],targetBarycentricLocation[1],targetBarycentricLocation[2]);
    };

//Code copy-pastes a lot of the displace particle routine above...eventually refactor more nicely
void triangulatedMeshSpace::transportParticleAndVelocity(meshPosition &pos, vector3 &velocityVector, vector3 &displacementVector)
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
    //dual check to ensure we are actually in the tangent plane with both direction of motion and velocity vector
    if (abs(currentSourceNormal*displacementVector) > THRESHOLD)
        {
        displacementVector -= (currentSourceNormal*displacementVector)*currentSourceNormal;
        }
    if (abs(currentSourceNormal*velocityVector) > THRESHOLD)
        {
        velocityVector -= (currentSourceNormal*velocityVector)*currentSourceNormal;
        }
    point3 target = sourcePoint + displacementVector;
    vector3 currentMove = vector3(sourcePoint, target);
    halfedgeIndex lastUsedHalfedge(-1);

int iter= 0;
    while(continueShifting)
        {
iter+=1;
        //get the current barycentric coordinates of the target
        targetBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0],vertexPositions[1],vertexPositions[2],target);
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
	
	if(uninvolvedVertex.size()== 2)
            {
            if (surface.is_border(vertexIndex(involvedVertex[0])))
                {
                targetBarycentricLocation = intersectionPoint;
                clampToThreshold(targetBarycentricLocation); //clamp avoids any near-edge confusion at later steps
                continueShifting = false;
                continue;
                }
            toIntersection = vector3(sourcePoint, PMP::construct_point(std::make_pair(currentSourceFace,intersectionPoint), surface)); //always assumed to be from source to intersected vertex
            std::pair<faceIndex, vector3> targetHeading = throughVertex(involvedVertex[0], toIntersection, currentSourceFace);
            provisionalTargetFace = targetHeading.first;
            newHeading = targetHeading.second; //this heading comes out normalized
            }
        if(uninvolvedVertex.size() != 1)
            ERRORERROR("a barycentric coordinate of the target is negative, but neither 1 nor 2 intersections were found. Apparently some debugging is needed!");
        /*
        Assume at this point that only one intersection point was found.
        intersectionPoint contains the barycentric coordinates of the intersection point
        on a current face. Identify the next face from the relevant halfEdge, update current
        move to be from the edge intersection to the original target, and rotate it around
        the edge according to the angle of the face normals.
        */
        point3 edgeIntersectionPoint = globalSMSP->point(currentSourceFace,intersectionPoint);
        sourcePoint = edgeIntersectionPoint;

        //update the move vector to intersection->target
        currentMove = vector3(edgeIntersectionPoint, target);
        vector3 targetNormal;
        halfedgeIndex intersectedEdge;

	if (uninvolvedVertex.size() == 1)
            {
            //if we're going through an edge, standard rotation around the axis formed by the normals is fine
            intersectedEdge = surface.halfedge(involvedVertex[0],involvedVertex[1]);
            //as in displaceParitcle, particle stops completely if it hits a boundary, for now
            if (surface.is_border(edgeIndex(intersectedEdge)))
                {
                targetBarycentricLocation = intersectionPoint;
                clampToThreshold(targetBarycentricLocation);
                continueShifting = false;
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
            axisVector /= vectorMagnitude(axisVector); //this could use normalize if we wanted to
            std::vector<point3> axis = {sourcePoint, sourcePoint+axisVector};
            target = rotateAboutAxis(target, axis, angle);
            displacementVector = vector3(sourcePoint, target); //move gets updated here
            //we can only update last used halfedge when crossing an edge -- maybe ask daniel if this breaks for going through vertices
            lastUsedHalfedge = intersectedEdge;
            //rotate the velocity vector in the same way
            point3 velocityVectorTarget = sourcePoint + velocityVector;
            velocityVectorTarget=rotateAboutAxis(velocityVectorTarget,axis,angle);
            velocityVector =  vector3(sourcePoint,velocityVectorTarget);
            }

	//two conditions are exclusive
        if (uninvolvedVertex.size() == 2)
            {
            //through a vertex, we use the straightness criterion of poltier et al 2006
            double displacementLength = vectorMagnitude(displacementVector);
            displacementVector = newHeading*displacementLength;
            targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
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
	    //finally, set intersected edge to either of the edges containing the involved vertex so it doesn't fail to update
	    intersectedEdge = surface.halfedge(involvedVertex[0],vertexList[uninvolvedVertex[0]]); 
            }

        //update source face index info and new bary coords in  the new face.
        currentSourceFace = provisionalTargetFace;
        getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
        sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0],vertexPositions[1],vertexPositions[2],edgeIntersectionPoint);
        currentSourceNormal = targetNormal;
        lastUsedHalfedge = intersectedEdge; //we probably need to do this using the vertex we went through
        
	if(iter ==maximumShiftEdgeCrossings ) ERRORERROR("Error: shifted across an extremely  large number of faces... is this an error, or an accidental call with an extremely large displacement vector? ");
        };

    pos.faceIndex = currentSourceFace;
    pos.x = point3(targetBarycentricLocation[0],targetBarycentricLocation[1],targetBarycentricLocation[2]);
    };

