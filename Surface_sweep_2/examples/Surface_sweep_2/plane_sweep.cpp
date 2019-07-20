//! \file examples/Arrangement_on_surface_2/plane_sweep.cpp
// Computing intersection points among curves using the surface-sweep alg.

#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Curve_2                                       Segment_2;

int main()
{
  // Construct the input segments.
  Segment_2 segments[] = {Segment_2 (Point_2 (1, 5), Point_2 (8, 5)),
                          Segment_2 (Point_2 (1, 1), Point_2 (8, 8)),
                          Segment_2 (Point_2 (3, 1), Point_2 (3, 8)),
                          Segment_2 (Point_2 (8, 5), Point_2 (8, 8))};

  // Compute all intersection points.
  std::list<Point_2> pts;

  CGAL::compute_intersection_points(segments, segments + 4,
                                    std::back_inserter(pts));

  // Print the result.
  std::cout << "Found " << pts.size() << " intersection points: " << std::endl;
  std::copy(pts.begin(), pts.end(),
            std::ostream_iterator<Point_2>(std::cout, "\n"));

  // Compute the non-intersecting sub-segments induced by the input segments.
  std::list<Segment_2> sub_segs;

  CGAL::compute_subcurves(segments, segments + 4, std::back_inserter(sub_segs));

  std::cout << "Found " << sub_segs.size()
            << " interior-disjoint sub-segments." << std::endl;

  CGAL_assertion(CGAL::do_curves_intersect (segments, segments + 4));

  return 0;
}
