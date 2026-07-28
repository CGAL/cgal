#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Snap_rounding_2.h>
#include <CGAL/Double_grid_snap_rounding_traits_2.h>

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

  // Compute the snapped subsegments and check if they do intersect
  std::vector< Polyline_2 > out;
  CGAL::snap_rounding_2(segs, out, CGAL::parameters::geom_traits(CGAL::Double_grid_snap_rounding_traits_2<Kernel>()).output_unique_segments(true));
  std::vector< Segment_2 > out_segs;
  out_segs.reserve(out.size());
  for(const Polyline_2& pl: out)
    out_segs.emplace_back(pl[0], pl[1]);
  std::cout << "Does the output intersect: " << CGAL::Surface_sweep_2::do_intersect(out_segs.begin(), out_segs.end(), false) << std::endl;

  return 0;
}
