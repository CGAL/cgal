// Implementation: Testprogram for 3-dimensional Segment Trees
// A three dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.


#include <CGAL/Segment_tree_k.h>
#include "include/Tree_Traits.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <list>

typedef CGAL::Segment_tree_3<CGAL::Tree_traits_3> Segment_tree_3_type;

int main()
{
  typedef CGAL::Tree_traits_3::Interval Interval;
  typedef CGAL::Tree_traits_3::Key Key;
  // definition of the two-dimensional segment tree
  std::list<Interval> InputList, OutputList, N;

  // insertion of the tree elements into the sequence container
  InputList.push_back(Interval(Key(1, 3, 5), Key(2,5,7)));
  InputList.push_back(Interval(Key(2,6.7,5), Key(4,6.9, 8)));
  InputList.push_back(Interval(Key(2,4.55, 8), Key(5, 7.88, 10)));
  InputList.push_back(Interval(Key(2, 4.66, 5), Key(6, 8.99, 8)));

  // creation of the segment tree
  std::list<Interval>::iterator first = InputList.begin();
  std::list<Interval>::iterator last = InputList.end();

  Segment_tree_3_type Segment_tree_3(first,last);

  // perform a window query
  Interval a(Key(3,4.88, 6), Key(6, 8.999, 9));
  Segment_tree_3.window_query(a,std::back_inserter(OutputList));

  // output of the querey elements on stdout
  std::list<Interval>::iterator j = OutputList.begin();
  std::cerr << "\n window_query (3,4.88, 6),(6, 8.999, 9)\n";
  while(j!=OutputList.end())
  {
    std::cerr << (*j).first.key_1 << "," << (*j).first.key_2 << ", "
	 << (*j).first.key_3 << "-" << (*j).second.key_1 << ","
	 << (*j).second.key_2 << ", " << (*j).second.key_3 << std::endl;
    j++;
  }
  Interval b(Key(2,6.8,9),Key(3,7,10));
  Segment_tree_3.enclosing_query(b,std::back_inserter(N));
  j = N.begin();
  std::cerr << "\n enclosing_query(2,6.8,9),(3,7,10)\n";
  while(j!=N.end())
  {
    std::cerr << (*j).first.key_1 << "," << (*j).first.key_2 << ", "
	 << (*j).first.key_3 << "-" << (*j).second.key_1 << ","
	 << (*j).second.key_2 << ", " << (*j).second.key_3 << std::endl;
    j++;
  }
  if(Segment_tree_3.segment_tree_3->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0;
}
