#ifndef ARR_CIRCULAR_H
#define ARR_CIRCULAR_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel;
typedef Kernel::FT                                         Number_type;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>          Traits;
typedef Traits::CoordNT                                    CoordNT;
typedef Traits::Point_2                                    Point;
typedef Traits::Curve_2                                    Curve;
typedef Traits::X_monotone_curve_2                         X_monotone_curve;
typedef Traits::Rational_point_2                           Rational_point;
typedef Traits::Rational_segment_2                         Segment;
typedef Traits::Rational_circle_2                          Circle;
typedef CGAL::Arrangement_2<Traits>                        Arrangement;

#endif
