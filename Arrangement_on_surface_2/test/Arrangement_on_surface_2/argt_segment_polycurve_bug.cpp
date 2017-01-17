// Constructing an arrangement of polylines.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>
/*
  Define the Arrangement traits class to be used. You can either use some user
  defined kernel and Segment_traits_2 or the defaults.
 */
// Instantiate the traits class using a user-defined kernel
// and Segment_traits_2.
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>     Geom_traits_2;
// Identical instantiation can be achieved using the default Kernel:
// typedef CGAL::Arr_polyline_traits_2<>                    Geom_traits_2;
typedef Geom_traits_2::Point_2                            Point_2;
typedef Geom_traits_2::Segment_2                          Segment_2;
typedef Geom_traits_2::Curve_2                            Polyline_2;
typedef Geom_traits_2::X_monotone_curve_2                 XM_2;
typedef CGAL::Arrangement_2<Geom_traits_2>                Arrangement_2;
int main()
{
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);
  Geom_traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();

  std::vector<Point_2> points1;
  points1.push_back(Point_2(0, 0));
  points1.push_back(Point_2(320, 0));
  std::vector<Point_2> points2;
  points2.push_back(Point_2(16, 0));
  points2.push_back(Point_2(64, 32));
  points2.push_back(Point_2(96, 0));
  points2.push_back(Point_2(144, 0));
  points2.push_back(Point_2(192, 32));
  points2.push_back(Point_2(224, 0));
  points2.push_back(Point_2(288, 0));
  std::vector<Point_2> points3;
  points3.push_back(Point_2(0, -16));
  points3.push_back(Point_2(64, 48));
  points3.push_back(Point_2(192, -32));
  std::vector<Polyline_2> pls(3);
  pls[0] = polyline_construct(points1.begin(), points1.end());
  pls[1] = polyline_construct(points2.begin(), points2.end());
  pls[2] = polyline_construct(points3.begin(), points3.end());
  CGAL::insert(arr, pls.begin(), pls.end());
  return 0;
}
