//! \file examples/Arrangement_on_surface_2/polyline_construction_tests.cpp
//  Testing the various possible constructions of polylines.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

typedef CGAL::Quotient<CGAL::MP_Float>                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef Segment_traits_2::Curve_2                       Segement_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Polyline_2;
typedef Traits_2::X_monotone_curve_2                    X_polyline_2;

int main ()
{
  // TODO: Add test with isolated points.
  // TODO: Test with long chains of sub-polylines

  Traits_2 traits;
  Traits_2::Construct_curve_2 polyline_const =
    traits.construct_curve_2_object();
  Traits_2::Construct_x_monotone_curve_2 x_polyline_const =
    traits.construct_x_monotone_curve_2_object();

  // Two polylines that will be constructed over and over again
  Polyline_2 poly;
  X_polyline_2 x_poly;


  // Construction from two points
  Point_2 p1 = Point_2(0,0);
  Point_2 p2 = Point_2(1,1);
  // poly = polyline_const(p1,p2);
  // x_poly = x_polyline_const(p1,p2);

  // Construction of x-mono from range of points
  std::cout << "Starting a construction of x-monotone from a range of point"
            << std::endl;
  std::vector<Point_2> pts;
  pts.push_back(Point_2(2,-11));
  pts.push_back(Point_2(1,0));
  // pts.push_back(Point_2(0,-1));
  x_poly = x_polyline_const(pts.begin(),pts.end());
  std::cout << "The constructed polyline is:" << std::endl;
  std::cout << x_poly << std::endl;

  // // Simple construction test of x-mono polyline from a range of segments.
  // std::vector<Segement_2> segs1;
  // segs1.push_back(Segement_2(Point_2(0,0),Point_2(1,1)));
  // segs1.push_back(Segement_2(Point_2(1,1),Point_2(5,2)));
  // X_polyline_2 xpoly1 = x_polyline_const(segs1.begin(),segs1.end());
  return 0;
}
