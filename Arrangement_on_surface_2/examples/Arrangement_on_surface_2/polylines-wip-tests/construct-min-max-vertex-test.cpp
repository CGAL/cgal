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

  std::vector<Point_2> pts;
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(1,1));
  pts.push_back(Point_2(3,-1));
  pts.push_back(Point_2(3.1,1));

  X_Polyline_2 xpoly = xpolyline_construct(pts.begin(),pts.end());

  std::cout << "For the xpoly: " << xpoly << "\n"
            << "The MIN vertex is: "
            << traits.construct_min_vertex_2_object()(xpoly) << "\n"
            << "The MAX vertex is: "
            << traits.construct_max_vertex_2_object()(xpoly) << "\n";

  xpoly = xpolyline_construct(pts.rbegin(),pts.rend());
  std::cout << "For the reversed xpoly: " << xpoly << "\n"
            << "The MIN vertex is: "
            << traits.construct_min_vertex_2_object()(xpoly) << "\n"
            << "The MAX vertex is: "
            << traits.construct_max_vertex_2_object()(xpoly) << "\n";

  pts.clear();
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(0,1));
  pts.push_back(Point_2(0,2));
  pts.push_back(Point_2(0,3));
  xpoly = xpolyline_construct(pts.begin(),pts.end());
  std::cout << "For the vertical xpoly: " << xpoly << "\n"
            << "The MIN vertex is: "
            << traits.construct_min_vertex_2_object()(xpoly) << "\n"
            << "The MAX vertex is: "
            << traits.construct_max_vertex_2_object()(xpoly) << "\n";

  xpoly = xpolyline_construct(pts.rbegin(),pts.rend());
  std::cout << "For the reversed vertical xpoly: " << xpoly << "\n"
            << "The MIN vertex is: "
            << traits.construct_min_vertex_2_object()(xpoly) << "\n"
            << "The MAX vertex is: "
            << traits.construct_max_vertex_2_object()(xpoly) << "\n";

  return 0;
}
