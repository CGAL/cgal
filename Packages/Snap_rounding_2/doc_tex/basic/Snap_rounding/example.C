#include "Isr_2.h"

typedef leda_real Number_Type;

typedef CGAL::Segment_2<CGAL::Cartesian<Number_Type> > Segment_2;
typedef CGAL::Point_2<CGAL::Cartesian<Number_Type> > Point_2;

int main()
{
  list<Segment_2> seg_list;

  seg_list.push_back(Segment_2(Point_2(0,0),Point_2(10,10)));
  seg_list.push_back(Segment_2(Point_2(0,10),Point_2(10,0)));
  seg_list.push_back(Segment_2(Point_2(3,0),Point_2(3,10)));
  seg_list.push_back(Segment_2(Point_2(7,0),Point_2(7,10)));

  ISR<Number_Type> i(seg_list.begin(),seg_list.end(),1.0);

  int counter = 0;
  for(list<list<Point_2> >::const_iterator iter1 = i.begin();iter1 != i.end();++iter1) {
    cout << "Polyline number " << ++counter << ":\n";
    for(list<Point_2>::const_iterator iter2 = iter1->begin();iter2 != iter1->end();++iter2)
      cout << "    (" << iter2->x().to_double() << ":" << iter2->y().to_double() << ")\n";
  }

  return(0);
}


