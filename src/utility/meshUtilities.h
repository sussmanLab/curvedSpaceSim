#ifndef meshUtilities_h
#define meshUtilities_h

#include "cgalIncludesAndTypedefs.h"
#include "pointDataType.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/boost/graph/iterator.h>

typedef CGAL::Surface_mesh<point3>                                      triangleMesh;
typedef triangleMesh::Face_index                                        faceIndex;
typedef triangleMesh::Halfedge_index                                    halfedgeIndex;
typedef triangleMesh::Vertex_index                                      vertexIndex;

namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                pmpBarycentricCoordinates;
typedef PMP::Face_location<triangleMesh, FT>                            pmpFaceLocation;
typedef CGAL::Face_around_target_circulator<triangleMesh>               faceCirculator;
typedef typename boost::graph_traits<triangleMesh>::vertex_descriptor   vertexDescriptor;
typedef typename boost::graph_traits<triangleMesh>::face_descriptor     faceDescriptor;

typedef boost::graph_traits<triangleMesh>                               graphTraits;
typedef graphTraits::vertex_iterator                                    vertexIterator;
typedef graphTraits::face_iterator                                      faceIterator;

typedef CGAL::Surface_mesh_shortest_path_traits<K, triangleMesh>        Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits>                        surfaceMeshShortestPath;
typedef typename surfaceMeshShortestPath::Barycentric_coordinates       smspBarycentricCoordinates;
typedef typename surfaceMeshShortestPath::Face_location                 smspFaceLocation;
typedef surfaceMeshShortestPath::Shortest_path_result                   shortestPathResult;

typedef CGAL::AABB_face_graph_triangle_primitive<triangleMesh>          AABB_face_graph_primitive;
typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>            AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                         AABB_tree;

smspFaceLocation meshPositionToFaceLocation(const meshPosition &p);
void getVertexPositionsFromFace(triangleMesh &mesh, faceIndex i, std::vector<point3> &result);
void getVertexIndicesFromFace(triangleMesh &mesh, faceIndex i, std::vector<point3> &result);

#endif
