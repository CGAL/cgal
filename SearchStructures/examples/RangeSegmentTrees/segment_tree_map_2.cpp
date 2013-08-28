// Implementation: Testprogram for 2-dimensional Segment Trees
// A two dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.

#include <iostream>
#include <utility>
#include <vector>
#include <iterator>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_tree_k.h>
#include <CGAL/Range_segment_tree_traits.h>



typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Segment_tree_map_traits_2<K, char> Traits;
typedef CGAL::Segment_tree_2<Traits > Segment_tree_2_type;

int main()
{
  typedef Traits::Interval Interval;
  typedef Traits::Pure_interval Pure_interval;
  typedef Traits::Key Key;
  // definition of the two-dimensional segment tree
  std::list<Interval> InputList, OutputList, N;



  // insertion of the tree elements into the sequence container
  InputList.push_back(Interval(Pure_interval(Key(1,5), Key(2,7)),'a'));
  InputList.push_back(Interval(Pure_interval(Key(2,7), Key(3,8)),'b'));
  InputList.push_back(Interval(Pure_interval(Key(6,9), Key(9,13)),'c'));
  InputList.push_back(Interval(Pure_interval(Key(1,3), Key(3,9)),'d'));

  // creation of the segment tree
  std::list<Interval>::iterator first = InputList.begin();
  std::list<Interval>::iterator last = InputList.end();

  Segment_tree_2_type Segment_tree_2(first,last);

  // perform a window query
  Interval a=Interval(Pure_interval(Key(3,6), Key(7,12)),'e');
  Segment_tree_2.window_query(a,std::back_inserter(OutputList));

  // output of the querey elements on stdout
  std::list<Interval>::iterator j = OutputList.begin();
  std::cerr << "\n window_query (3,6),(7,12)\n";
  while(j!=OutputList.end())
  {
    std::cerr << (*j).first.first.x() << "-" << (*j).first.second.x() << " "
	 << (*j).first.first.y() << "-" << (*j).first.second.y() << std::endl;
    j++;
  }



  Interval b=Interval(Pure_interval(Key(6,10),Key(7,11)), 'f');
  Segment_tree_2.enclosing_query(b,std::back_inserter(N));
  j = N.begin();
  std::cerr << "\n enclosing_query (6,10),(7,11)\n";
  while(j!=N.end())
  {
    std::cerr << (*j).first.first.x() << "-" << (*j).first.second.x() << " "
	 << (*j).first.first.y() << "-" << (*j).first.second.y() << std::endl;
    j++;
  }
  if(Segment_tree_2.segment_tree_2->is_valid())
    std::cerr << "Tree  is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0;
}
