#ifndef CGAL_INCLUDES_TYPEDEFS_H
#define CGAL_INCLUDES_TYPEDEFS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

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
#endif
