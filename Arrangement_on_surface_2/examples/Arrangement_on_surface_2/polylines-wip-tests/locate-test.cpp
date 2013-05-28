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

const bool TO_RIGHT = true;
const bool TO_LEFT = false;

void test_locate(std::vector<Point_2> pts, Point_2 p)
{
  Traits_2 traits;
  Traits_2::Construct_x_monotone_curve_2 xpolyline_construct =
    traits.construct_x_monotone_curve_2_object();

  std::cout << "** Test the point (" << p << ")\n";

  std::cout << "*** Direct Case\n";
  X_Polyline_2 xpoly = xpolyline_construct(pts.begin(),pts.end());
  std::cout << "Consider: " << xpoly << "\n";

  int res = traits.locate(xpoly,p);
  std::cout << "locate(xpoly,p) = " << res << "\n";
  res = traits.locate_side(xpoly,p,TO_RIGHT);
  std::cout << "locate_side(xpoly,p,TO_RIGHT) = " << res << "\n";
  res = traits.locate_side(xpoly,p,TO_LEFT);
  std::cout << "locate_side(xpoly,p,TO_LEFT) = " << res << "\n";

  std::cout << "*** Inverse case\n";
  xpoly = xpolyline_construct(pts.rbegin(),pts.rend());
  std::cout << "Consider: " << xpoly << "\n";
  res = traits.locate(xpoly,p);
  std::cout << "locate(cv,p) = " << res << "\n";
  res = traits.locate_side(xpoly,p,TO_RIGHT);
  std::cout << "locate_side(xpoly,p,TO_RIGHT) = " << res << "\n";
  res = traits.locate_side(xpoly,p,TO_LEFT);
  std::cout << "locate_side(xpoly,p,TO_LEFT) = " << res << "\n";
}

int main ()
{

  Traits_2 traits;
  Traits_2::Construct_x_monotone_curve_2 xpolyline_construct =
    traits.construct_x_monotone_curve_2_object();

  std::cout << "* Vertical xpoly with one segment\n";

  std::vector<Point_2> pts;
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(0,1));

  // Contains
  test_locate(pts,Point_2(0,0));
  test_locate(pts,Point_2(0,1));
  test_locate(pts,Point_2(0,0.5));
  // Above/below
  test_locate(pts,Point_2(0,2));
  test_locate(pts,Point_2(0,-1));
  // Left/right
  test_locate(pts,Point_2(0.2,0));
  test_locate(pts,Point_2(-0.2,0));
  test_locate(pts,Point_2(0.2,1));
  test_locate(pts,Point_2(-0.2,1));
  test_locate(pts,Point_2(0.2,0.5));
  test_locate(pts,Point_2(-0.2,0.5));
  test_locate(pts,Point_2(0.2,-1));
  test_locate(pts,Point_2(-0.2,-1));
  test_locate(pts,Point_2(0.2,2));
  test_locate(pts,Point_2(-0.2,2));

  std::cout << "* Vertical xpoly with two segment\n";
  pts.push_back(Point_2(0,2));
  test_locate(pts,Point_2(0,0));
  test_locate(pts,Point_2(0,1));
  test_locate(pts,Point_2(0,2));
  test_locate(pts,Point_2(0,0.5));
  test_locate(pts,Point_2(0,1.4));
  test_locate(pts,Point_2(0.5,0));
  test_locate(pts,Point_2(-0.5,0));
  test_locate(pts,Point_2(0.5,1));
  test_locate(pts,Point_2(-0.5,1));
  test_locate(pts,Point_2(0.5,2));
  test_locate(pts,Point_2(-0.5,2));
  test_locate(pts,Point_2(0.5,1.2));
  test_locate(pts,Point_2(-0.5,1.2));
  test_locate(pts,Point_2(0.5,-1));
  test_locate(pts,Point_2(-0.5,2.1));

  std::cout << "* Vertical xpoly with more segments\n";
  pts.push_back(Point_2(0,3));
  pts.push_back(Point_2(0,4));
  pts.push_back(Point_2(0,5));
  pts.push_back(Point_2(0,6));

  test_locate(pts,Point_2(0,0));
  test_locate(pts,Point_2(0,6));
  test_locate(pts,Point_2(0,1));
  test_locate(pts,Point_2(0,5));
  test_locate(pts,Point_2(0,3));
  test_locate(pts,Point_2(0,0.5));
  test_locate(pts,Point_2(0,5.6));
  test_locate(pts,Point_2(0,1.6));
  test_locate(pts,Point_2(0,3.6));
  test_locate(pts,Point_2(0,6));
  test_locate(pts,Point_2(0.5,0));
  test_locate(pts,Point_2(-0.5,0));
  test_locate(pts,Point_2(0.5,1));
  test_locate(pts,Point_2(-0.5,1));
  test_locate(pts,Point_2(0.5,6));
  test_locate(pts,Point_2(-0.5,6));
  test_locate(pts,Point_2(0.5,4.1));
  test_locate(pts,Point_2(-0.5,4.1));
  test_locate(pts,Point_2(0.5,-1));
  test_locate(pts,Point_2(-0.5,-1));
  test_locate(pts,Point_2(0.5,7));
  test_locate(pts,Point_2(-0.5,7));

  std::cout << "* Test non-vertical xpoly with one segment\n";
  pts.clear();
  pts.push_back(Point_2(0,0));
  pts.push_back(Point_2(1,1));
  test_locate(pts,Point_2(-1,-2));
  test_locate(pts,Point_2(-1,-1));
  test_locate(pts,Point_2(-1,2));
  test_locate(pts,Point_2(0,-2));
  test_locate(pts,Point_2(0,0));
  test_locate(pts,Point_2(0,2));
  test_locate(pts,Point_2(0.5,-1));
  test_locate(pts,Point_2(0.5,0.5));
  test_locate(pts,Point_2(0.5,1));
  test_locate(pts,Point_2(1,-1));
  test_locate(pts,Point_2(1,1));
  test_locate(pts,Point_2(1,2));
  test_locate(pts,Point_2(2,-1));
  test_locate(pts,Point_2(2,2));
  test_locate(pts,Point_2(2,2.2));

  std::cout << "* Test non-vertical xpoly with two segments\n";
  pts.push_back(Point_2(2,2));
  test_locate(pts,Point_2(-1,-2));
  test_locate(pts,Point_2(-1,-1));
  test_locate(pts,Point_2(-1,2));
  test_locate(pts,Point_2(0,-2));
  test_locate(pts,Point_2(0,0));
  test_locate(pts,Point_2(0,2));
  test_locate(pts,Point_2(0.5,-1));
  test_locate(pts,Point_2(0.5,0.5));
  test_locate(pts,Point_2(0.5,1));
  test_locate(pts,Point_2(1,-1));
  test_locate(pts,Point_2(1,1));
  test_locate(pts,Point_2(1,2));
  test_locate(pts,Point_2(2,-1));
  test_locate(pts,Point_2(2,2));
  test_locate(pts,Point_2(2,2.2));
  test_locate(pts,Point_2(3,-2.2));
  test_locate(pts,Point_2(3,3));
  test_locate(pts,Point_2(3,6));

  std::cout << "* Test non-vertical xpoly with more segments\n";
  pts.push_back(Point_2(3,3));
  pts.push_back(Point_2(4,4));
  pts.push_back(Point_2(5,5));
  test_locate(pts,Point_2(-1,-2));
  test_locate(pts,Point_2(-1,-1));
  test_locate(pts,Point_2(-1,2));
  test_locate(pts,Point_2(0,-2));
  test_locate(pts,Point_2(0,0));
  test_locate(pts,Point_2(0,2));
  test_locate(pts,Point_2(0.5,-1));
  test_locate(pts,Point_2(0.5,0.5));
  test_locate(pts,Point_2(0.5,1));
  test_locate(pts,Point_2(1,-1));
  test_locate(pts,Point_2(1,1));
  test_locate(pts,Point_2(1,2));
  test_locate(pts,Point_2(5,4));
  test_locate(pts,Point_2(5,5));
  test_locate(pts,Point_2(5,6));
  test_locate(pts,Point_2(6,5));
  test_locate(pts,Point_2(6,6));
  test_locate(pts,Point_2(6,7));
  //TODO: Verify and complete the test of locte(cv,p) for a multi-segments cv


  return 0;
}
