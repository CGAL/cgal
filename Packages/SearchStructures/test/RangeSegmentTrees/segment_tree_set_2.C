
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_tree_k.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <iostream>
#include <utility>
#include <vector>
#include <iterator>
#include <list>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Range_segment_tree_set_traits_2<K> Traits;
typedef CGAL::Segment_tree_2<Traits > Segment_tree_2_type;

typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef std::pair<Point_2, Point_2> Interval;
int main()
{

  std::list<Iso_rectangle_2> rectangles;
  std::copy(std::istream_iterator<Iso_rectangle_2>(std::cin),
	    std::istream_iterator<Iso_rectangle_2>(),
	    std::back_inserter(rectangles));

  std::list<Interval> intervals;
  for(std::list<Iso_rectangle_2>::iterator it = rectangles.begin();
      it != rectangles.end();
	++it){
    intervals.push_back(Interval(it->vertex(0),it->vertex(2)));
  }

  std::list<Interval>  output, N;

  // creation of the segment tree
  std::list<Interval>::iterator first = intervals.begin();
  std::list<Interval>::iterator last = intervals.end();

  Segment_tree_2_type Segment_tree_2(first, last);

  // perform a window query
  Interval a=Interval(Point_2(3,6), Point_2(7,12));
  Segment_tree_2.window_query(a,std::back_inserter(output));

  // output of the querey elements on stdout
  std::list<Interval>::iterator j = output.begin();
  std::cerr << "\n window_query (3,6), (7,12)\n";
  while(j!=output.end())
  {
    std::cerr << (*j).first.x() << "-" << (*j).second.x() << " " 
	 << (*j).first.y() << "-" << (*j).second.y() << std::endl; 
    j++;
  }
  std::cerr << "\n enclosing_query (6,10),(7,11) \n";
  Interval b=Interval(Point_2(6,10),Point_2(7,11));
  Segment_tree_2.enclosing_query(b,std::back_inserter(N));
  j = N.begin();
  std::cerr << "\n enclosing_query (6,10),(7,11) \n";
  while(j!=N.end())
  {
    std::cerr << (*j).first.x() << "-" << (*j).second.x() << " " 
	 << (*j).first.y() << "-" << (*j).second.y() << std::endl; 
    j++;
  }
  if(Segment_tree_2.segment_tree_2->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0;
}


