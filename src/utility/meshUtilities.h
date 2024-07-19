#ifndef meshUtilities_h
#define meshUtilities_h

#include "cgalIncludesAndTypedefs.h"
#include "pointDataType.h"
#include "functionUtilities.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <memory>

typedef CGAL::Surface_mesh<point3>                                      triangleMesh;
typedef triangleMesh::Face_index                                        faceIndex;
typedef triangleMesh::Halfedge_index                                    halfedgeIndex;
typedef triangleMesh::Vertex_index                                      vertexIndex;
typedef triangleMesh::Edge_index                                        edgeIndex;

namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                pmpBarycentricCoordinates;
typedef PMP::Face_location<triangleMesh, FT>                            pmpFaceLocation;
typedef CGAL::Face_around_target_circulator<triangleMesh>               faceCirculator;
typedef typename boost::graph_traits<triangleMesh>::vertex_descriptor   vertexDescriptor;
typedef typename boost::graph_traits<triangleMesh>::face_descriptor     faceDescriptor;
typedef CGAL::Vertex_around_target_circulator<triangleMesh>             vertexCirculator;

typedef boost::graph_traits<triangleMesh>                               graphTraits;
typedef graphTraits::vertex_iterator                                    vertexIterator;
typedef graphTraits::face_iterator                                      faceIterator;

typedef CGAL::Surface_mesh_shortest_path_traits<K, triangleMesh>        Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits>                        surfaceMeshShortestPath;
typedef typename surfaceMeshShortestPath::Barycentric_coordinates       smspBarycentricCoordinates;
typedef typename surfaceMeshShortestPath::Face_location                 smspFaceLocation;
typedef surfaceMeshShortestPath::Shortest_path_result                   shortestPathResult;

typedef CGAL::AABB_face_graph_triangle_primitive<triangleMesh>          AABB_face_graph_primitive;
typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>                 AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                         AABB_tree;

typedef CGAL::Face_filtered_graph<triangleMesh>                         filteredGraph;

smspFaceLocation meshPositionToFaceLocation(const meshPosition &p);
//!REQUIRE result has size 3
void getVertexPositionsFromFace(triangleMesh &mesh, faceIndex i, std::vector<point3> &result);
//!REQUIRE result has size 3
void getVertexIndicesFromFace(triangleMesh &mesh, const faceIndex &i, std::vector<vertexIndex> &result);

double meanEdgeLength(triangleMesh mesh, bool verbose = false);
double triangleArea(point3 v1, point3 v2, point3 v3);
double meanTriangleArea(triangleMesh mesh);
double totalArea(triangleMesh mesh); 

void clampToThreshold(pmpBarycentricCoordinates &baryPoint);

//!return true if the two lines which pass through the given endpoints intersect between the specified points on both lines. fill in the barycentric location of the intersection point
bool intersectionOfLinesInBarycentricCoordinates(pmpBarycentricCoordinates line1Start, pmpBarycentricCoordinates line1End, pmpBarycentricCoordinates line2Start, pmpBarycentricCoordinates line2End, pmpBarycentricCoordinates &intersectionPoint);


//specialize the intersectionOfLinesInBarycentricCoordinates function to the case where you care about the intersection of "line2" with a "line1" which goes from one vertex to another (i.e., has barycentric coordinates which are a permutation of (1,0,0))
bool intersectionBarycentricLinesV1V2(pmpBarycentricCoordinates line2Start, pmpBarycentricCoordinates line2End, pmpBarycentricCoordinates &intersectionPoint);
bool intersectionBarycentricLinesV2V3(pmpBarycentricCoordinates line2Start, pmpBarycentricCoordinates line2End, pmpBarycentricCoordinates &intersectionPoint);
bool intersectionBarycentricLinesV3V1(pmpBarycentricCoordinates line2Start, pmpBarycentricCoordinates line2End, pmpBarycentricCoordinates &intersectionPoint);

//Go through each of the vAvB edge intersection functions
bool findTriangleEdgeIntersectionInformation(pmpBarycentricCoordinates sourceBarycentricLocation, pmpBarycentricCoordinates targetBarycentricLocation, pmpBarycentricCoordinates &intersectionPoint, std::vector<vertexIndex> vertexList, halfedgeIndex  previousHalfEdge, triangleMesh &surface, std::vector<vertexIndex> &involvedVertex,std::vector<int> &uninvolvedVertex);

void convertBarycentricCoordinates(triangleMesh &mesh1, triangleMesh &mesh2, std::unordered_map<faceIndex,int> &faceMap, smspFaceLocation &locationToConvert);

void computePathDistanceAndTangents(surfaceMeshShortestPath *smsp, smspFaceLocation &targetPoint, double &distance, vector3 &startPathTangent, vector3 &endPathTangent);

void printPoint(point3 a);
void printBary(smspBarycentricCoordinates a);

#endif
