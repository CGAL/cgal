#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/hot_pixel_snap_rounding_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                          Point_2;
typedef Kernel::Segment_2                        Segment_2;

typedef std::vector<Point_2>                                    Polyline_2;

int main()
{
  std::vector<Segment_2> segs;

  segs.push_back(Segment_2(Point_2(0, 0), Point_2(10, 10)));
  segs.push_back(Segment_2(Point_2(0, 10), Point_2(10, 0)));
  segs.push_back(Segment_2(Point_2(3, 0), Point_2(3, 10)));
  segs.push_back(Segment_2(Point_2(7, 0), Point_2(7, 10)));

  // Compute the snapped subsegments and check if they do intersect
  std::vector< Polyline_2 > out;
  CGAL::hot_pixel_snap_rounding_2(segs, out, CGAL::parameters::do_iterative_snap_rounding(false));
  std::vector< Segment_2 > out_segs;
  out_segs.reserve(out.size());
  for(const Polyline_2& pl: out)
    out_segs.emplace_back(pl[0], pl[1]);
  std::cout << "Does the output intersect: " << CGAL::do_curves_intersect(out_segs.begin(), out_segs.end()) << std::endl;

  return 0;
}
