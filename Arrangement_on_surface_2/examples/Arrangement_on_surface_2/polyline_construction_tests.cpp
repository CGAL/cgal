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
typedef Segment_traits_2::Curve_2                       Segment_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Polyline_2;
typedef Traits_2::X_monotone_curve_2                    X_polyline_2;

int main ()
{
  Traits_2 traits;
  Traits_2::Construct_curve_2 polyline_const =
    traits.construct_curve_2_object();
  Traits_2::Construct_x_monotone_curve_2 x_polyline_const =
    traits.construct_x_monotone_curve_2_object();

  // Two polylines that will be constructed over and over again
  Polyline_2 poly;
  X_polyline_2 x_poly;

  // Containers of points and segments
  std::vector<Point_2> pts;

  std::cout << "Starting a construction of a polyline from two points..."
            << std::endl;
  Point_2 p1 = Point_2(1,1);
  Point_2 p2 = Point_2(0,0);
  poly = polyline_const(p1,p2);
  x_poly = x_polyline_const(p1,p2);
  std::cout << "The constructed polyline is:" << std::endl;
  std::cout << "non-x-mono: " << poly << std::endl;
  std::cout << "x-mono    : " << x_poly << std::endl;
  std::cout << "----====----"<< std::endl;

  std::cout << "Starting a construction of a polyline from a range of point..."
            << std::endl;
  pts.push_back(Point_2(2,-11));
  pts.push_back(Point_2(1,0));
  pts.push_back(Point_2(0,-1));
  pts.push_back(Point_2(0,1));
  pts.push_back(Point_2(-1,-1));
  pts.push_back(Point_2(0,1));
  poly = polyline_const(pts.begin(),pts.end());
  pts.clear();
  std::cout << "The constructed polyline is:" << std::endl;
  std::cout << poly << std::endl;
  std::cout << "----====----"<< std::endl;

  // Construction of polyline from range of points
  // std::cout << "Starting a construction of polyline from a range of point"
  //           << std::endl;
  // pts.clear();
  // pts.push_back(Point_2(0,0));
  // pts.push_back(Point_2(1,1));
  // pts.push_back(Point_2(2,0));
  // pts.push_back(Point_2(2,-3));
  // pts.push_back(Point_2(-5,0));
  // pts.push_back(Point_2(20,0));
  // poly = polyline_const(pts.begin(),pts.end());
  // std::cout << "The constructed polyline is:" << std::endl;
  // std::cout << poly << std::endl;
  // std::cout << "----====----"<< std::endl;

  // Construction from a single segment
  // std::cout << "Starting a construction of polyline from a single segment."
  //           << std::endl;
  // Segment_2 seg = Segment_2(Point_2(0,0),Point_2(1,1));
  // poly = polyline_const(seg);
  // std::cout << "The constructed polyline is:" << std::endl;
  // std::cout << poly << std::endl;
  // std::cout << "----====----"<< std::endl;

  return 0;
}
