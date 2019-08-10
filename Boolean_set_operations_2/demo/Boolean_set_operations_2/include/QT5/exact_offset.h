#ifndef EXACT_OFFSET_H
#define EXACT_OFFSET_H

#include <fstream>
#include <CGAL/config.h>
#include <boost/timer.hpp>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>         Traits;
typedef CGAL::General_polygon_set_2<Traits>               Polygon_set_2;
typedef Traits::Polygon_2                                 Polygon_2_12;
typedef Traits::Polygon_with_holes_2                      Polygon_with_holes_2_12;

template <typename Polygon_2> Polygon_with_holes_2_12 app_offsetting(const Polygon_2_12 & P, const typename Kernel::FT & r)
{
	Polygon_with_holes_2_12 offset = CGAL::approximated_offset_2(P,r,0.00001);
	return offset;
}

#endif