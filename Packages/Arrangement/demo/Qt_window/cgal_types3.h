#ifndef CGAL_TYPES_HEADER
#define CGAL_TYPES_HEADER

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
//#include <CGAL/Largest_empty_iso_rectangle_2.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>  Coord_type;
typedef CGAL::Cartesian<Coord_type> Rep;

typedef Rep::Point_2             Point;
typedef Rep::Segment_2           Segment;
//typedef Rep::Iso_rectangle_2     Iso_rectangle_2;


typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Seg_traits;
typedef CGAL::Arr_polyline_traits_2<Seg_traits>         Traits;

//typedef Traits::Point_2                                 Point;
typedef Seg_traits::Curve_2                             Curve;
typedef Seg_traits::X_monotone_curve_2                  X_curve;

typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arr_base_node<Traits::Curve, Traits::X_curve> Base_node;
typedef CGAL::Arrangement_2<Dcel, Traits, Base_node>    Arr;

typedef CGAL::Polygon_traits_2<Rep>         PT;
typedef std::vector<Point>                  Container;
typedef CGAL::Polygon_2<PT, Container>      Cgal_Polygon;

#endif
