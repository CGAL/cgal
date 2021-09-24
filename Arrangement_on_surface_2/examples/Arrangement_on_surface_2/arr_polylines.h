#ifndef ARR_POLYLINES_H
#define ARR_POLYLINES_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT                                        Number_type;

typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>       Traits;
typedef Traits::Point_2                                   Point;
typedef Traits::Segment_2                                 Segment;
typedef Traits::Curve_2                                   My_polyline;
typedef CGAL::Arrangement_2<Traits>                       Arrangement;

#endif
