#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>           Number_type;
typedef CGAL::Cartesian<Number_type>             Kernel;
typedef CGAL::Snap_rounding_traits_2<Kernel>     Traits;
typedef Kernel::Segment_2                        Segment_2;
typedef Kernel::Point_2                          Point_2;
typedef std::vector<Segment_2>                   Segment_range_2;
typedef std::vector<Point_2>                     Polyline_2;
typedef std::vector<Polyline_2>                  Polyline_range_2;

int main()
{
  Segment_range_2 segs;
  Polyline_range_2 output;

  segs.push_back(Segment_2(Point_2(0, 0), Point_2(10, 10)));
  segs.push_back(Segment_2(Point_2(0, 10), Point_2(10, 0)));
  segs.push_back(Segment_2(Point_2(3, 0), Point_2(3, 10)));
  segs.push_back(Segment_2(Point_2(7, 0), Point_2(7, 10)));

  // Generate an iterated snap-rounding representation, where the centers of
  // the hot pixels bear their original coordinates, using 5 kd trees:
  // CGAL::snap_rounding_2<Traits,Segment_list_2::const_iterator,Polyline_list_2>
  //   (seg_list.begin(), seg_list.end(), output_list, 1.0, true, false, 5);

  CGAL::hot_pixel_snap_rounding_2(segs, std::back_inserter(output));
  int counter = 0;
  for (const Polyline_2 &pl: output) {
    std::cout << "Polyline number " << ++counter << ":\n";
    Polyline_2::const_iterator iter2;
    for (const Point_2 &p: pl)
      std::cout << "    (" << p.x() << ":" << p.y() << ")\n";
  }

  return(0);
}
