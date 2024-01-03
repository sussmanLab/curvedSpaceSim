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
    };

void triangulatedMeshSpace::distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent)
    {
    int nTargets = p2.size();
    distances.resize(nTargets);
    startPathTangent.resize(nTargets);
    endPathTangent.resize(nTargets);

    surfaceMeshShortestPath pathFinder(surface);
    AABB_tree tree;
    pathFinder.build_aabb_tree(tree);

    //shortestPaths needs barycentric coordinates... for now this entails a conversion step. Potentially motivates a switch so that the meshPosition data type is always this faceLocation structure?
    faceLocation sourcePoint = pathFinder.locate<AABB_face_graph_traits>(p1.x,tree);
    pathFinder.add_source_point(sourcePoint.first,sourcePoint.second);
  
    pathFinder.build_sequence_tree();

    for(int ii = 0; ii < nTargets; ++ii)
        {
        //pathPoints holds the sequence of intersection points between the shortest path and the meshed surface (edges, vertices, etc)
        std::vector<point3> pathPoints; 
        faceLocation targetPoint = pathFinder.locate<AABB_face_graph_traits>(p2[ii].x,tree);
        shortestPathResult geodesic = pathFinder.shortest_path_points_to_source_points(targetPoint.first, targetPoint.second,  std::back_inserter(pathPoints));
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
    };

void triangulatedMeshSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector)
    {
    UNWRITTENCODE("all mesh space stuff unwritten");
    };
