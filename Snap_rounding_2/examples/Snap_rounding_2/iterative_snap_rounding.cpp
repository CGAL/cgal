#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/hot_pixel_snap_rounding_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Segment_2                        Segment_2;
typedef Kernel::Point_2                          Point_2;
typedef std::vector<Point_2>                     Polyline_2;
typedef std::vector<Polyline_2>                  PolylineRange;

int main()
{
  std::vector<Segment_2> segs;

  segs.push_back(Segment_2(Point_2(0, 0), Point_2(10, 10)));
  segs.push_back(Segment_2(Point_2(0, 10), Point_2(10, 0)));
  segs.push_back(Segment_2(Point_2(3, 0), Point_2(3, 10)));
  segs.push_back(Segment_2(Point_2(7, 0), Point_2(7, 10)));

  // Generate an iterated snap-rounding representation, where the centers of
  // the hot pixels are integers, using one kd tree, which is the default:
  std::vector<Polyline_2> out;
  CGAL::hot_pixel_snap_rounding_2(segs, std::back_inserter(out));

  int counter = 0;
  for(const Polyline_2 &pl: out) {
    std::cout << "Polyline number " << ++counter << ":" << std::endl;
    for(const Point_2 &p: pl)
      std::cout << p << std::endl;
  }

  return 0;
}
