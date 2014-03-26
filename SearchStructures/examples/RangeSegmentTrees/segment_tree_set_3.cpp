// Implementation: Testprogram for 3-dimensional Segment Trees
// A three dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.

#include <iostream>
#include <utility>
#include <iterator>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_tree_k.h>
#include <CGAL/Range_segment_tree_traits.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Range_segment_tree_set_traits_3<K> Traits;
typedef CGAL::Segment_tree_3<Traits > Segment_tree_3_type;

int main()
{
  typedef Traits::Interval Interval;
  typedef Traits::Key Key;
  // definition of the two-dimensional segment tree
  std::list<Interval> InputList, OutputList, N;

  // insertion of the tree elements into the sequence container
  InputList.push_back(Interval(Key(1,5,7), Key(2,7,9)));
  InputList.push_back(Interval(Key(2,7,6), Key(3,8,9)));
  InputList.push_back(Interval(Key(6,9,5), Key(9,13,8)));
  InputList.push_back(Interval(Key(1,3,4), Key(3,9,8)));

  // creation of the segment tree
  std::list<Interval>::iterator first = InputList.begin();
  std::list<Interval>::iterator last = InputList.end();

  Segment_tree_3_type Segment_tree_3(first,last);

  // perform a window query
  Interval a(Key(3,6,5), Key(7,12,8));
  Segment_tree_3.window_query(a,std::back_inserter(OutputList));

  // output of the querey elements on stdout
  std::list<Interval>::iterator j = OutputList.begin();
  std::cerr << "\n window_query (3,6,5),(7,12,8) \n";
  while(j!=OutputList.end())
  {
    std::cerr << (*j).first.x() << "," << (*j).first.y() << "," << (*j).first.z()
	 <<", " << (*j).second.x() << "," << (*j).second.y() << ","
	 << (*j).second.z() << std::endl;
    j++;
  }
  Interval b(Key(6,10,7),Key(7,11,8));
  Segment_tree_3.enclosing_query(b,std::back_inserter(N));
  j = N.begin();
  std::cerr << "\n enclosing_query (6,10,7), (7,11,8)\n";
  while(j!=N.end())
  {
    std::cerr << (*j).first.x() << "," << (*j).first.y() << "," << (*j).first.z()
	 <<", " << (*j).second.x() << "," << (*j).second.y() << ","
	 << (*j).second.z() << std::endl;
    j++;
  }
  if(Segment_tree_3.segment_tree_3->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0;
}
