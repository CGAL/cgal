#ifndef BOPS_CIRCULAR_H
#define BOPS_CIRCULAR_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>         Traits;

typedef CGAL::General_polygon_set_2<Traits>               Polygon_set_2;
typedef Traits::Polygon_2                                 Polygon_2;
typedef Traits::Polygon_with_holes_2                      Polygon_with_holes_2;

#endif
