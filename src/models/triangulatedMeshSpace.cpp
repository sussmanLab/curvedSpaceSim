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
/*
    double distanceToTravel = vectorMagnitude(displacementVector);
    pmpFaceLocation sourceLocation = PMP::locate(pos.x,surface);
    point3 sourcePoint = PMP::construct_point(sourceLocation,surface);
    faceIndex currentSourceFace = sourceLocation.first;
    vector3 currentSourceNormal = PMP::compute_face_normal(currentSourceFace,surface);
    if (abs(currentSourceNormal*displacementVector) > 1e-18)
        {
        ERRORERROR("non-tangent displacement vector on a face");
        }

    point3 target = sourcePoint + displacementVector;
    std::vector<vertexIndex> vertexList;
//    std::vector<point3> vertexPositions = getVertexPositions(surface,sourceLocation.first);
*/
/*
  double travelLength = vectorMagnitude(move);

  //short pseudo-projection routine in case we're slightly off the face (so code is general to arbitrary starting point) 
  Face_location sourceLocation = PMP::locate(pos, mesh);
  Point_3 source_point = PMP::construct_point(sourceLocation, mesh); //xyz point representing current source

  Face_index currentSourceFace = sourceLocation.first;  
  Vector_3 currentSourceNormal = PMP::compute_face_normal(currentSourceFace,mesh);
 
  Point_3 target = source_point+move;

  std::vector<Vertex_index> vertexList;
  std::vector<Point_3> vertexPos = getVertexPositions(mesh, sourceLocation.first);

  //if the movement vector isn't quite in the tangent plane, put it there. 
  if (abs(currentSourceNormal*move) > 0) {
    target = project_to_face(vertexPos, pos+move);
  }
  Vector_3 current_move = Vector_3(source_point, target); 
  
  //if there is no intersection, avoid initalizing any of the intersection code
  std::array<double,3> targetBary = PMP::barycentric_coordinates(vertexPos[0],vertexPos[1],vertexPos[2],target);
  bool intersection = false;
  
  for (double bc: targetBary) {  
    if (bc < 0) {
      intersection = true; 
    }
  }

  if (intersection == false) {
     return source_point + current_move; 
  }

  //initializations if there *is* an intersection
  std::vector<Point_3> targetVertices;
  std::vector<Point_3> forRotation;
  Point_3 rotatedTarget;
  Face_index currentTargetFace;
  //Vector_3 currentTargetFaceNormal;
  double lengthToSharedElement;
  double rotationAngle;
  double overlap;

  while(intersection){
    //vertexList = getVertexPositions(mesh,currentSourceFace);
    vertexList = getVertexIndices(mesh,currentSourceFace);


    std::pair<Point_3, std::vector<Vertex_index>> intersection_info = find_intersection(mesh, currentSourceFace, source_point, source_point+current_move, vertexList);
    Point_3 intersection_point = intersection_info.first;
    std::vector<Vertex_index> intersected_elements = intersection_info.second;
    
    if (intersection_point == Point_3(1000,1000,1000)) {
      intersection = false;
      break;
    }

    Vector_3 vector_to_intersection = Vector_3(source_point, intersection_point);
    double lengthToSharedElement = vectorMagnitude(vector_to_intersection); //how far we've traveled
    current_move = reduceVector(current_move, lengthToSharedElement);         //decrease move size by length to intersected vertex/edge --
                                                                            //effectively the step where we "walk" to that intersection
    
    target = intersection_point+current_move; //storage of where the move vector currently points for rotation later
    
    currentTargetFace = getTargetFace(intersected_elements, vector_to_intersection, currentSourceFace, mesh); //face we're about to walk into;
    
    source_point = intersection_point;//update source to be the most recent intersection point -- finish walking there

    Face_location newMoveLocation = rotateIntoNewFace(mesh, currentSourceFace, currentTargetFace, source_point, target);
    Point_3 rotatedTarget = PMP::construct_point(newMoveLocation,mesh);

    //check that we've rotated in the right direction via overlap
    current_move = Vector_3(source_point, rotatedTarget);//source is now intersection
    currentSourceFace = currentTargetFace;
  }
  //source_point+move is the location in the original face if there were no intersections, and it will 
  //be the location in the unfolded mesh if there were intersections (from an edge intersection to a spot
  //within a face)
  
  Point_3 fTarget = source_point+current_move;
  return fTarget;  
  */
UNWRITTENCODE("displace particle");
    };

void triangulatedMeshSpace::findIntersection(faceIndex sourceFace, point3 source, point3 target, std::vector<vertexIndex> &vertexIndices, point3 &intersectionPoint, std::vector<vertexIndex> &intersections)
    {
    UNWRITTENCODE("findIntersection unwritten");
    }

