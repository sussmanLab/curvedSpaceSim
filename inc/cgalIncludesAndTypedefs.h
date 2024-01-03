#ifndef CGAL_INCLUDES_TYPEDEFS_H
#define CGAL_INCLUDES_TYPEDEFS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
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

//types and names
typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::FT                                                           FT;
typedef K::Point_2                                                      Point_2;
typedef K::Vector_2                                                     Vector_2;
typedef K::Ray_2                                                        Ray_2;
typedef K::Segment_2                                                    Segment_2;
typedef K::Point_3                                                      point3;
typedef K::Vector_3                                                     vector3;
typedef K::Ray_3                                                        ray3;
typedef K::Segment_3                                                    segment3;
typedef K::Intersect_3                                                  intersect3;
typedef K::Triangle_3                                                   triangle3;
typedef CGAL::Surface_mesh<point3>                                      triangleMesh;
typedef triangleMesh::Face_index                                        faceIndex;
typedef triangleMesh::Halfedge_index                                    halfedgeIndex;
typedef triangleMesh::Vertex_index                                      vertexIndex;
typedef CGAL::Surface_mesh_shortest_path_traits<K, triangleMesh>        Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits>                        surfaceMeshShortestPath;


namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
//typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
//typedef PMP::Face_location<triangleMesh, FT>                            faceLocation;
typedef CGAL::Face_around_target_circulator<triangleMesh>               faceCirculator;
typedef typename boost::graph_traits<triangleMesh>::vertex_descriptor   vertexDescriptor;
typedef typename boost::graph_traits<triangleMesh>::face_descriptor     faceDescriptor;

typedef boost::graph_traits<triangleMesh>                               graphTraits;
typedef graphTraits::vertex_iterator                                    vertexIterator;
typedef graphTraits::face_iterator                                      faceIterator;

typedef typename surfaceMeshShortestPath::Barycentric_coordinates       barycentricCoordinates;
typedef typename surfaceMeshShortestPath::Face_location                 faceLocation;
typedef surfaceMeshShortestPath::Shortest_path_result                   shortestPathResult;

typedef CGAL::AABB_face_graph_triangle_primitive<triangleMesh>          AABB_face_graph_primitive;
typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>            AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                         AABB_tree;
#endif
