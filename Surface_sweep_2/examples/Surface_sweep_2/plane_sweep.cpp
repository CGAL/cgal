//! \file examples/Surface_sweep_2/plane_sweep.cpp
// Computing intersection points among curves using the surface-sweep alg.

#include <list>
#include <cassert>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Segment_2                                         Segment;
typedef Kernel::Point_2                                         Point;
using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
using My_polyline = Traits::Curve_2;

int main()
{
    Traits traits;
  // Construct the input segments.
  auto polyline_construct = traits.construct_curve_2_object();

  Point points1[5];
  points1[0] = Point(0, 0);
  points1[1] = Point(2, 4);
  points1[2] = Point(3, 0);
  points1[3] = Point(4, 4);
  points1[4] = Point(6, 0);
  auto pi1 = polyline_construct(&points1[0], &points1[5]);

  std::list<Point> points2;
  points2.push_back(Point(1, 3));
  points2.push_back(Point(0, 2));
  points2.push_back(Point(1, 0));
  points2.push_back(Point(2, 1));
  points2.push_back(Point(3, 0));
  points2.push_back(Point(4, 1));
  points2.push_back(Point(5, 0));
  points2.push_back(Point(6, 2));
  points2.push_back(Point(5, 3));
  points2.push_back(Point(4, 2));
  auto pi2 = polyline_construct(points2.begin(), points2.end());

  std::vector<Segment> segs;
  segs.push_back(Segment(Point(0, 2), Point(1, 2)));
  segs.push_back(Segment(Point(1, 2), Point(3, 6)));
  segs.push_back(Segment(Point(3, 6), Point(5, 2)));
  auto pi3 = polyline_construct(segs.begin(), segs.end());

  My_polyline vec[3]={pi1,  pi2, pi3};

  // Compute all intersection points.
  std::list<Point_2> pts;

  CGAL::compute_intersection_points(vec, vec + 3,
                                    std::back_inserter(pts));

  // Print the result.
  std::cout << "Found " << pts.size() << " intersection points: " << std::endl;
  std::copy(pts.begin(), pts.end(),
            std::ostream_iterator<Point_2>(std::cout, "\n"));

  // Compute the non-intersecting sub-segments induced by the input segments.
  std::list<My_polyline> sub_segs;

  CGAL::compute_subcurves(vec, vec + 3, std::back_inserter(sub_segs));

  std::cout << "Found " << sub_segs.size()
            << " interior-disjoint sub-segments." << std::endl;

//  assert(CGAL::do_curves_intersect (segments, segments + 4));

  return 0;
}
