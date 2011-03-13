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
typedef std::list<Segment_2>                     Segment_list_2;
typedef std::list<Point_2>                       Polyline_2;
typedef std::list<Polyline_2>                    Polyline_list_2;

int main()
{
  Segment_list_2 seg_list;
  Polyline_list_2 output_list;

  seg_list.push_back(Segment_2(Point_2(0, 0), Point_2(10, 10)));
  seg_list.push_back(Segment_2(Point_2(0, 10), Point_2(10, 0)));
  seg_list.push_back(Segment_2(Point_2(3, 0), Point_2(3, 10)));
  seg_list.push_back(Segment_2(Point_2(7, 0), Point_2(7, 10)));

  // Generate an iterated snap-rounding representation, where the centers of
  // the hot pixels bear their original coordinates, using 5 kd trees:
  CGAL::snap_rounding_2<Traits,Segment_list_2::const_iterator,Polyline_list_2>
    (seg_list.begin(), seg_list.end(), output_list, 1.0, true, false, 5);

  int counter = 0;
  Polyline_list_2::const_iterator iter1;
  for (iter1 = output_list.begin(); iter1 != output_list.end(); ++iter1) {
    std::cout << "Polyline number " << ++counter << ":\n";
    Polyline_2::const_iterator iter2;
    for (iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2)
      std::cout << "    (" << iter2->x() << ":" << iter2->y() << ")\n";
  }

  return(0);
}
