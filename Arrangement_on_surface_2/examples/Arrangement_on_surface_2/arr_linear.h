#ifndef ARR_LINEAR_H
#define ARR_LINEAR_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT                                        Number_type;

typedef CGAL::Arr_linear_traits_2<Kernel>                 Traits;
typedef Traits::Point_2                                   Point;
typedef Traits::Segment_2                                 Segment;
typedef Traits::Ray_2                                     Ray;
typedef Traits::Line_2                                    Line;
typedef Traits::X_monotone_curve_2                        X_monotone_curve;

typedef CGAL::Arrangement_2<Traits>                       Arrangement;
typedef Arrangement::Vertex_handle                        Vertex_handle;
typedef Arrangement::Halfedge_handle                      Halfedge_handle;
typedef Arrangement::Face_handle                          Face_handle;
typedef Arrangement::Vertex_const_handle                  Vertex_const_handle;
typedef Arrangement::Halfedge_const_handle                Halfedge_const_handle;
typedef Arrangement::Face_const_handle                    Face_const_handle;

#endif
