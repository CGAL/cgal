#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Gmpq.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Bso_segment_traits_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Bso_circle_segment_traits_2.h>

typedef CGAL::Gmpq                                    Coord_type;
//typedef double                                        Coord_type;
typedef CGAL::Cartesian<Coord_type>		                Kernel;

typedef Kernel::Segment_2						                  Segment;
typedef Kernel::Point_2						                    Point;
typedef Kernel::Circle_2                              Circle;
//typedef CGAL::Polygon_2<Kernel>			                  Polygon;
//typedef CGAL::General_polygon_with_holes_2<Polygon>   Polygon_with_holes_2;
//typedef Polygon_with_holes_2::Holes_const_iterator    Holes_const_iterator;


//typedef CGAL::Bso_segment_traits_2<Kernel>            Traits;
typedef CGAL::Bso_circle_segment_traits_2<Kernel>     Traits;
typedef Traits::Curve_2                               Curve;
typedef Traits::X_monotone_curve_2                    XCurve;
typedef Traits::Polygon_2                             Polygon;
typedef CGAL::General_polygon_with_holes_2<Polygon>   Polygon_with_holes_2;
typedef CGAL::General_polygon_set_2<Traits>           Polygon_set;
typedef Polygon_with_holes_2::Holes_const_iterator    Holes_const_iterator;



typedef CGAL::Iso_rectangle_2<Kernel>                 Iso_rectangle;
#endif
