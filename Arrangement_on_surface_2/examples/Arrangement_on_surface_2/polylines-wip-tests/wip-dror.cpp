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

#include "../arr_print.h"

typedef CGAL::Quotient<CGAL::MP_Float>                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Segment_traits_2::Curve_2                       Segment_2;
typedef Traits_2::Curve_2                               Polyline_2;
typedef Traits_2::X_monotone_curve_2                    X_polyline_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

int main ()
{

  Traits_2 traits;
  Traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();
  Traits_2::Construct_x_monotone_curve_2 x_polyline_construct =
    traits.construct_x_monotone_curve_2_object();

  Arrangement_2         arr;

  Point_2 p1;
  Point_2 p2;
  Point_2 p3;
  Point_2 p4;

  /* Test the functor Intersect_2 */
  std::list<Point_2>    pts;
  p1 = Point_2(0,0); pts.push_back(p1);
  p2 = Point_2(2,2); pts.push_back(p2);
  p3 = Point_2(4,0); pts.push_back(p3);
  X_polyline_2 poly1 = x_polyline_construct(pts.begin(), pts.end());
  std::cout << "poly1 = " << poly1 << std::endl;

  pts.clear();
  p1 = Point_2(0,1); pts.push_back(p1);
  p2 = Point_2(3,1); pts.push_back(p2);
  p3 = Point_2(5,1); pts.push_back(p3);
  X_polyline_2 poly2 = x_polyline_construct(pts.begin(), pts.end());
  std::cout << "poly2 = " << poly2 << std::endl;

  Traits_2::Intersect_2 int_2 = traits.intersect_2_object();
  std::vector<CGAL::Object> out;
  int_2(poly1,poly2,back_inserter(out));
  std::cout << "Number of intersection elements: " << out.size() << std::endl;
  for (
       std::vector<CGAL::Object>::iterator it = out.begin() ;
       it!=out.end() ; ++it)
    {
      typedef std::pair<Point_2,Traits_2::Multiplicity> pt_int;
      const pt_int p =
        CGAL::object_cast<pt_int>(*it);
      std::cout << "Intersection element: " << p.first << std::endl;
    }
  /* END: Test the functor Intersect_2 */


  /* Test Equal_2 of two x-polys */
  Traits_2::Equal_2 equal_2 = traits.equal_2_object();
  X_polyline_2 xcv1,xcv2;
  pts.clear();
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(0,0.2));
  pts.push_back(Point_2(0,0.3));
  pts.push_back(Point_2(0,0.4));
  pts.push_back(Point_2(0,0.5));
  pts.push_back(Point_2(0,0.75));
  xcv1 = x_polyline_construct(pts.begin(),pts.end());
  pts.clear();
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(0,0.5));
  pts.push_back(Point_2(0,0.75));
  xcv2 = x_polyline_construct(pts.begin(),pts.end());
  std::cout << "Testing Equal_2\n"
    "xcv1 is: " << xcv1 << "\n"
    "xcv2 is: " << xcv2 << "\n";
  std::cout << "Equal(xcv1,xcv2) = " << equal_2(xcv1,xcv2) << std::endl;
  std::cout << "Equal(xcv2,xcv1) = " << equal_2(xcv2,xcv1) << std::endl;
  /* END: Test Equal_2 of two x-polys */

  /* Test construction from range of segments */
  std::cout << "\nTest construction from range of segs\n";
  pts.clear();
  p1 = Point_2(0,0);
  p2 = Point_2(1,2);
  p3 = Point_2(3,0);
  p4 = Point_2(4,3);
  std::vector<Segment_2> segs;
  segs.push_back(Segment_2(p1,p2));
  segs.push_back(Segment_2(p2,p3));
  segs.push_back(Segment_2(p3,p4));
  Polyline_2 poly = polyline_construct(segs.begin(), segs.end());
  std::cout << "Construction of general poly from range of points:\n"
            << poly << std::endl;

  /* Test construction from two points */
  std::cout << "\nConstruction (x-mono) from two points:\n"
            << x_polyline_construct(Point_2(0,0),Point_2(1,0)) << std::endl;

  /* Test construction from a degenerated Kernel::segment */
  Kernel::Segment_2 s(Kernel::Point_2(-10,0),Kernel::Point_2(0,0));
  // Segment_2 s(Point_2(0,0), Point_2(0,0));
  std::cout << "\nConstruction (x-mono) from a single segment:\n"
            << x_polyline_construct(s) << std::endl;


  return 0;
}
