/*
 * Test the functors:
 * Construct_min_vertex_2 and Construct_max_vertex_2
 */
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
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Segment_traits_2::Curve_2                       Segment_2;
typedef Traits_2::Curve_2                               Polyline_2;
typedef Traits_2::X_monotone_curve_2                    X_Polyline_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

int main ()
{

  Traits_2 traits;
  Traits_2::Construct_x_monotone_curve_2 xpolyline_construct =
    traits.construct_x_monotone_curve_2_object();
  Traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();
  Traits_2::Push_back_2 push_back = traits.push_back_2_object();

  std::vector<Point_2> pts;
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(1,1));
  pts.push_back(Point_2(3,0));
  pts.push_back(Point_2(0,10));

  Polyline_2 poly = polyline_construct(pts.begin(),pts.end());
  std::cout << "Consider the polyline:\n" << poly << "\n";
  Point_2 p(5,5);
  std::cout << "Push the point (" << p << ") and obtain:\n";
  push_back(poly,p);
  std::cout << poly << "\n";

  return 0;
}
