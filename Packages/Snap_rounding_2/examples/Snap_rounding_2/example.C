#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include "../../include/CGAL/Snap_rounding_2.h"

#include <CGAL/leda_real.h>

typedef leda_rational                    Number_type;
typedef CGAL::Cartesian<Number_type>     Rep;
typedef CGAL::Snap_rounding_traits<Rep>  Sr_traits;
typedef CGAL::Snap_rounding_2<Sr_traits> Sr;
typedef Sr::Segment_2                    Segment_2;
typedef Sr::Point_2                      Point_2;
typedef Sr::Segments_container           Segments;

int main()
{
  Segments seg_list;

  seg_list.push_back(Segment_2(Point_2(0,0),Point_2(10,10)));
  seg_list.push_back(Segment_2(Point_2(0,10),Point_2(10,0)));
  seg_list.push_back(Segment_2(Point_2(3,0),Point_2(3,10)));
  seg_list.push_back(Segment_2(Point_2(7,0),Point_2(7,10)));

  Sr snap(seg_list.begin(),seg_list.end(),1.0);

  int counter = 0;
  for(Sr::Polyline_const_iterator iter1 = snap.polylines_begin();
      iter1 != snap.polylines_end();
      ++iter1) {
    std::cout << "Polyline number " << ++counter << ":\n";
    for(Sr::Point_const_iterator iter2 = iter1->begin();
        iter2 != iter1->end();
        ++iter2)
      std::cout << "    (" << iter2->x().to_double() << ":"
                << iter2->y().to_double() << ")\n";
  }

  return(0);
}
