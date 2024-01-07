#include "triangulatedMeshSpace.h"
/*! \file triangulatedMeshSpace.cpp */
#include <stdexcept>


void triangulatedMeshSpace::loadMeshFromFile(std::string filename, bool verbose)
    {
    if(verbose)
        printf("loading from file %s\n",filename.c_str());
    if(!CGAL::IO::read_polygon_mesh(filename, surface) || !CGAL::is_triangle_mesh(surface))
        {
        std::cerr << "Invalid input file." << std::endl;
        throw std::exception();
        };

    globalSMSP = make_shared<surfaceMeshShortestPath>(surface);
    AABB_tree globalTree;
    globalSMSP->build_aabb_tree(globalTree);
        
    };

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


void triangulatedMeshSpace::distanceWithSubmeshing(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent)
    {
    UNWRITTENCODE("distance with submeshing not implemented yet");
    //create submesh:



    //create new surfaceMeshShortestPath based on it, follow the logic of the main distance routine. refactor code so it doesn't repeat
    };


/*
As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex)
*/
void triangulatedMeshSpace::distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent)
    {
    if(submeshingActivated)
        {
        distanceWithSubmeshing(p1,p2,distances,startPathTangent,endPathTangent);
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
        //pathPoints holds the sequence of intersection points between the shortest path and the meshed surface (edges, vertices, etc)
        std::vector<point3> pathPoints; 
        shortestPathResult geodesic = globalSMSP->shortest_path_points_to_source_points(targetPoint.first, targetPoint.second,  std::back_inserter(pathPoints));
        distances[ii] = std::get<0>(geodesic);
        //Note that the path goes from the target to source, so if we want to know path tangent at the source for force calculation, we must use the *end* of points[]
        int pathSize = pathPoints.size();
        startPathTangent[ii] = vector3(pathPoints[pathSize-2],pathPoints[pathSize-1]);
        endPathTangent[ii] = vector3(pathPoints[0],pathPoints[1]);
        //normalize path tangents
        double normalization = sqrt(startPathTangent[ii].squared_length());
        startPathTangent[ii] /= normalization;
        normalization = sqrt(endPathTangent[ii].squared_length());
        endPathTangent[ii] /= normalization;
        };
    globalSMSP->remove_all_source_points();
    };

/*
As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex)
*/
void triangulatedMeshSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector)
    {
    double distanceToTravel = vectorMagnitude(displacementVector);
    pmpFaceLocation sourceLocation = meshPositionToFaceLocation(pos);
    point3 sourcePoint = PMP::construct_point(sourceLocation,surface);
    faceIndex currentSourceFace = sourceLocation.first;

    std::vector<vertexIndex> vertexList;
    std::vector<point3> vertexPositions;

    getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
    pmpBarycentricCoordinates sourceBarycentricLocation = sourceLocation.second;
    pmpBarycentricCoordinates targetBarycentricLocation;
    bool continueShifting = true;
    vector3 currentSourceNormal = PMP::compute_face_normal(currentSourceFace,surface);
    if (abs(currentSourceNormal*displacementVector) > THRESHOLD)
        {
        ERRORERROR("non-tangent displacement vector on a face");
        //target = project_to_face(vertexPos, sourcePoint+displacementVector);
        }
    point3 target = sourcePoint + displacementVector;
    vector3 currentMove = vector3(sourcePoint, target); 

    while(continueShifting)
        {
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
        pmpBarycentricCoordinatesintersectionPoint;
        std::vector<int> uninvolvedVertex;
        std::vector<vertexIndex> involvedVertex;

        findTriangleEdgeIntersectionInformation(sourceBarycentricLocation,targetBarycentricLocation,intersectionPoint, vertexList, involvedVertex,uninvolvedVertex);
        if(uninvolvedVertex.size()== 2)
            {UNWRITTENCODE("line goes through a vertex...write this routine");}
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
printf("{");
printPoint(sourcePoint);printf(",");
printPoint(target);printf(",");
printPoint(vertexPositions[0]);printf(",");
printPoint(vertexPositions[1]);printf(",");
printPoint(vertexPositions[2]);printf(",");
printPoint(edgeIntersectionPoint); printf("},");



        sourcePoint = edgeIntersectionPoint;
        //vector3 vectorToIntersection = vector3(sourcePoint,edgeIntersectionPoint);
        //double distanceToIntersectionPoint = sqrt(vectorToIntersection.squared_length());

        //update the move vector to intersection->target
        currentMove = vector3(edgeIntersectionPoint, target);

        //identify the next faceIndex from the shared intersected edge
        halfedgeIndex intersectedEdge = surface.halfedge(involvedVertex[0],involvedVertex[1]);
        faceIndex provisionalTargetFace = surface.face(intersectedEdge);
        if (provisionalTargetFace == currentSourceFace) 
            provisionalTargetFace= surface.face(surface.opposite(intersectedEdge));

        //find the new target after rotating current move into the tangent plane ofthe new face
        vector3 targetNormal = PMP::compute_face_normal(provisionalTargetFace,surface);
        double normalDotProduct = currentSourceNormal*targetNormal;
        double angle = acos(normalDotProduct);
        vector3 axisVector = CGAL::cross_product(currentSourceNormal,targetNormal);
        axisVector /= vectorMagnitude(axisVector);
        std::vector<point3> axis = {sourcePoint, sourcePoint+axisVector};
        target = rotateAboutAxis(target, axis,angle);
        displacementVector = vector3(sourcePoint, target);

        //update source face index info and new bary coords in  the new face.
        currentSourceFace = provisionalTargetFace;
        getVertexPositionsFromFace(surface,currentSourceFace, vertexPositions);
        sourceBarycentricLocation = PMP::barycentric_coordinates(vertexPositions[0],vertexPositions[1],vertexPositions[2],sourceBarycentricLocation);
        currentSourceNormal = targetNormal;
        };

    pos.faceIndex = sourceLocation.first;
    pos.x = point3(targetBarycentricLocation[0],targetBarycentricLocation[1],targetBarycentricLocation[2]);
    };
