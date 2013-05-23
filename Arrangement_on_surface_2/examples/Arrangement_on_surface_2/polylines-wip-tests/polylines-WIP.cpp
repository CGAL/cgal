//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

#include "arr_print.h"

typedef CGAL::Quotient<CGAL::MP_Float>                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Segment_traits_2::Curve_2                       Segment_2;
typedef Traits_2::Curve_2                               Polyline_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

int main ()
{

  Traits_2 traits;
  Traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();

  Arrangement_2         arr;

  std::list<Point_2>    pts;
  Point_2 p1(0,2);
  Point_2 p2(-1,-1);
  Point_2 p3(1,-1);
  pts.push_back(p1);
  pts.push_back(p2);
  pts.push_back(p3);
  pts.push_back(p1);
  Polyline_2 poly1 = polyline_construct(pts.begin(), pts.end());

  std::vector<Segment_2> segs;
  Point_2 q1(-1,1);
  Point_2 q2(0,-2);
  Point_2 q3(1,1);

  segs.push_back(Segment_2(q1,q3));
  segs.push_back(Segment_2(q1,q2));
  segs.push_back(Segment_2(q2,q3));
  Polyline_2 poly2 = polyline_construct(segs.begin(), segs.end());

  // TODO: Make sure to construct all possible of polylines. That is, from
  //       ranges of points, ranges of segments, ranges with only one segment
  //       or only two points, etc.

  insert (arr, poly1);
  insert (arr, poly2);

  // print_arrangement (arr);

  std::cout << "\nSummary:\n"
            << "V = " << arr.number_of_vertices() << "\n"
            << "E = " << arr.number_of_edges() << "\n"
            << "F = " << arr.number_of_faces() << std::endl;

  return 0;
}
