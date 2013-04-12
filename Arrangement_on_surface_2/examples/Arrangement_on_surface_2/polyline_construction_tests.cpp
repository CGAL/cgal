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

  std::cout << "--\n"
    "This example construct various polylines using the various"
    "construction functors that are provided in the Arr_polyline_traits_2\n";

  Traits_2 traits;
  Traits_2::Construct_curve_2 polyline_const =
    traits.construct_curve_2_object();
  Traits_2::Construct_x_monotone_curve_2 x_polyline_const =
    traits.construct_x_monotone_curve_2_object();

  // Two polylines that will be constructed over and over again
  Polyline_2 poly;
  X_polyline_2 x_poly;

  std::cout << "----==Construction from two points==----"<< std::endl;
  Point_2 p1 = Point_2(0,0);
  Point_2 p2 = Point_2(1,1);

  poly = polyline_const(p1,p2);
  std::cout << "Polyline is: " << poly << std::endl;

  x_poly = x_polyline_const(p1,p2);
  std::cout << "x-mono polyline is: " << x_poly << std::endl;
  x_poly = x_polyline_const(p2,p1);
  std::cout << "x-mono polyline is the same even if the points' order is "
            << "reversed: " << x_poly << std::endl;

  std::cout << "\n----==Construction from a single segment==----"<< std::endl;
  std::cout << "Polyline is: " << polyline_const(Segment_2(p1,p2))
            << std::endl;
  std::cout << "Reverting the points' order yields: "
            << polyline_const(Segment_2(p2,p1)) << std::endl;
  std::cout << "For x-monotone construction the order doesn't matter:\n"
            << x_polyline_const(Segment_2(p1,p2)) << "\n"
            << x_polyline_const(Segment_2(p1,p2)) << std::endl;

  std::cout << "\n----==Construction from a range of point==----"<< std::endl;
  std::cout << "* general polyline:" << std::endl;
  // Containers of points and segments
  std::vector<Point_2> pts;
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(1,0));
  pts.push_back(Point_2(3,3));
  pts.push_back(Point_2(1,2));
  pts.push_back(Point_2(4,1));
  pts.push_back(Point_2(4,3));
  std::cout << "The constructed polyline is:" << std::endl;
  std::cout <<  polyline_const(pts.begin(),pts.end()) << std::endl;

  std::cout << "* X-monotone case:" << std::endl;
  pts.clear();
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(1,1));
  pts.push_back(Point_2(3,-10));
  std::cout << "The constructed x-monotone polyline is independent of the "
            << "input's order:" << std::endl;
  std::cout <<  x_polyline_const(pts.begin(),pts.end()) << std::endl;
  std::cout <<  x_polyline_const(pts.rbegin(),pts.rend()) << std::endl;

  std::cout << "* Same goes, for example, for vertical polyline" << std::endl;
  pts.clear();
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(0,1));
  pts.push_back(Point_2(0,10));
  std::cout << "The constructed x-monotone polyline is independent of the "
            << "input's order:" << std::endl;
  std::cout <<  x_polyline_const(pts.begin(),pts.end()) << std::endl;
  std::cout <<  x_polyline_const(pts.rbegin(),pts.rend()) << std::endl;

  std::cout << "\n----==Construction from a range of segments==----\n"
    "* general polyline:" << std::endl;
  std::vector<Segment_2> segs;
  Point_2 q1(Point_2(0,0));
  Point_2 q2(Point_2(1,1));
  Point_2 q3(Point_2(2,0));
  Point_2 q4(Point_2(3,1));
  //TODO: Discussion: The orientation of each segment in the range seems not to
  //      play a role. This is good, isn't it? The polyline is defined by its
  //      segments, regardless of their order. Verify that the constructions
  //      (i.e. general and x-mono) are consistent with the decision.
  segs.push_back(Segment_2(q2,q1));
  segs.push_back(Segment_2(q2,q3));
  segs.push_back(Segment_2(q3,q4));
  std::cout << "Constructed polyline is:\n"
            << polyline_const(segs.begin(),segs.end()) << std::endl;
  std::cout << "Constructed x-monotone polyline is:\n"
            << x_polyline_const(segs.begin(),segs.end()) << std::endl;
  std::cout << "And in reversed order, constructed x-monotone polyline is:\n"
            << x_polyline_const(segs.rbegin(),segs.rend()) << std::endl;

  return 0;
}
