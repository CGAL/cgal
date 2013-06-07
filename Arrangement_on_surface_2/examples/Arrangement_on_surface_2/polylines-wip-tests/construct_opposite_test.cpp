/*
 * Test the functor:
 * Compare_endpoints_xy_2
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

void test_pts(std::vector<Point_2> pts) {
  Traits_2 traits;
  Traits_2::Construct_x_monotone_curve_2 xpolyline_construct =
    traits.construct_x_monotone_curve_2_object();
  X_Polyline_2 xpoly = xpolyline_construct(pts.begin(),pts.end());
  std::cout << "The opposite xpoly of xpoly: " << xpoly << "\n"
            << "is:\n"
            << traits.construct_opposite_2_object()(xpoly) << "\n";
  xpoly = xpolyline_construct(pts.rbegin(),pts.rend());
  std::cout << "The opposite xpoly of xpoly: " << xpoly << "\n"
            << "is:\n"
            << traits.construct_opposite_2_object()(xpoly) << "\n";
}

int main ()
{


  std::vector<Point_2> pts;
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(1,1));
  pts.push_back(Point_2(3,-1));
  pts.push_back(Point_2(3.1,1));
  test_pts(pts);

  pts.clear();
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(0,1));
  pts.push_back(Point_2(0,2));
  pts.push_back(Point_2(0,3));
  test_pts(pts);

  pts.clear();
  test_pts(pts);

  return 0;
}
