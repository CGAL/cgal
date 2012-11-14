#ifndef ARR_INEXACT_CONSTRUCTION_SEGMENTS_H
#define ARR_INEXACT_CONSTRUCTION_SEGMENTS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT                                          Number_type;

typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>      Traits;
typedef Traits::Point_2                                     Point;
typedef Traits::X_monotone_curve_2                          Segment;

typedef CGAL::Arrangement_2<Traits>                         Arrangement;
typedef Arrangement::Vertex_handle                          Vertex_handle;
typedef Arrangement::Halfedge_handle                        Halfedge_handle;
typedef Arrangement::Face_handle                            Face_handle;

#endif
