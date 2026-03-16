#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Float_snap_rounding_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Segment_2                                       Segment_2;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::FT                                              FT;

typedef std::vector<Point_2>                                    Polyline_2;

int main()
{
  // Example with a non-trivial rounding
  std::vector< Segment_2 > segs;
  FT e(std::pow(2, -60));

  segs.emplace_back(Point_2(1-e, 1), Point_2(-1-e, -1+2*e));
  segs.emplace_back(Point_2(e/2, e/2), Point_2(1, -1));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(2, 0));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(-2+e, -4));
  segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));

  // Compute the snapped subcurves and check if they do intersect
  std::vector< Segment_2 > out;
  CGAL::compute_snapped_subcurves_2(segs.begin(), segs.end(), std::back_inserter(out));
  std::cout << "Does the output intersect: " << CGAL::do_curves_intersect(out.begin(), out.end()) << std::endl;
  std::cout << "Size of the output: " << out.size() << std::endl;

  // Compute the polylines corresponding to the output rounded on single precision floats
  std::vector< Polyline_2 > polylines;
  CGAL::float_snap_rounding_2(segs.begin(), segs.end(), std::back_inserter(polylines), CGAL::parameters::geom_traits(CGAL::Float_snap_rounding_traits_2<Kernel>()));

  return 0;
}
