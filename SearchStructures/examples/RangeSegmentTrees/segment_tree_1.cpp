// Implementation: Testprogram for 1-dimensional Segment Trees
// A one dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.

#include "include/Tree_Traits.h"
#include <CGAL/Segment_tree_k.h>
#include <iostream>
#include <iterator>
#include <list>
#include <utility>

typedef CGAL::Segment_tree_1<CGAL::Tree_traits_1> SSegment_tree_1_type;

int main()
{
  typedef CGAL::Tree_traits_1::Interval Interval;
  // definition of the one-dimensional segment tree
  std::list<Interval> InputList, OutputList, N;

  // insertion of the tree elements into the sequence container
  InputList.push_back(Interval(64, 81));
  InputList.push_back(Interval(465, 499));
  InputList.push_back(Interval(288, 379));
  InputList.push_back(Interval(314, 375));

  // creation of the segment tree
  typedef std::list<Interval>::iterator l_iterator;
  l_iterator first = InputList.begin();
  l_iterator last = InputList.end();

  SSegment_tree_1_type segment_tree_1(first,last);

  // perform a window query
  Interval a(45,200);
  segment_tree_1.window_query(a,std::back_inserter(OutputList));

  // output of the querey elements on stdout
  l_iterator j = OutputList.begin();
  std::cerr << "\n window_query (45,200)\n";
  while(j!=OutputList.end())
  {
    std::cerr << (*j).first << "-" << (*j).second << std::endl;
    j++;
  }
  Interval b(320, 370);
  segment_tree_1.enclosing_query(b,std::back_inserter(N));
  j = N.begin();
  std::cerr << "\n enclosing_query (320, 370)\n";
  while(j!=N.end())
  {
    std::cerr << (*j).first << "-" << (*j).second << std::endl;
    j++;
  }
  if(segment_tree_1.segment_tree_1->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0;
}
