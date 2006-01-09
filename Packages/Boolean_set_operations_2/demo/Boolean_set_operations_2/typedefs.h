#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Gmpq.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>

//typedef CGAL::Quotient<CGAL::MP_Float>                Coord_type;
typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >  Coord_type;
//typedef double                                        Coord_type;
typedef CGAL::Cartesian<Coord_type>		                Kernel;

typedef Kernel::Segment_2						                  Segment;
typedef Kernel::Point_2						                    Point;
typedef Kernel::Circle_2                              Circle;

typedef CGAL::Gps_circle_segment_traits_2<Kernel>     Traits;
typedef Traits::Curve_2                               Curve;
typedef Traits::X_monotone_curve_2                    XCurve;
typedef Traits::Polygon_2                             Polygon;
typedef CGAL::General_polygon_with_holes_2<Polygon>   Polygon_with_holes;
typedef CGAL::General_polygon_set_2<Traits>           Polygon_set;
typedef Polygon_with_holes::Holes_const_iterator      Holes_const_iterator;

typedef CGAL::Iso_rectangle_2<Kernel>                 Iso_rectangle;
#endif
