// Implementation: Testprogram for 2-dimensional Segment Trees
// A two dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.


#include <CGAL/Segment_tree_k.h>
#include "include/Tree_Traits.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <list>

typedef CGAL::Segment_tree_2<CGAL::Tree_traits_2> Segment_tree_2_type;

int main()
{

  typedef CGAL::Tree_traits_2::Interval Interval;
  typedef CGAL::Tree_traits_2::Key Key;
  // definition of the two-dimensional segment tree
  std::list<Interval> InputList, OutputList, N;


  // insertion of the tree elements into the sequence container
  // InputList.push_back(Interval(Key(64, 139), Key(81.3, 921.3)));
  // InputList.push_back(Interval(Key(465, 504), Key(499.3, 829.0)));
  // InputList.push_back(Interval(Key(288, 875), Key(379.5, 982.7)));
  // InputList.push_back(Interval(Key(314, 465), Key(375.1, 711.5)));
  InputList.push_back(Interval(Key(64, 139), Key(81, 921.3)));
  InputList.push_back(Interval(Key(465, 504), Key(499, 829.0)));
  InputList.push_back(Interval(Key(288, 875), Key(379, 982.7)));
  InputList.push_back(Interval(Key(314, 465), Key(375, 711.5)));

  // creation of the segment tree
  std::list<Interval>::iterator first = InputList.begin();
  std::list<Interval>::iterator last = InputList.end();

  Segment_tree_2_type Segment_tree_2(first,last);

  // perform a window query
  // Interval a(Key(45, 500), Key(200.0, 675.1));
  Interval a(Key(45, 500), Key(200, 675.1));
  Segment_tree_2.window_query(a,std::back_inserter(OutputList));

  // output of the querey elements on stdout
  std::list<Interval>::iterator j = OutputList.begin();
  std::cerr << "\n window_query (45, 500), (200.0, 675.1)\n";
  while(j!=OutputList.end())
  {
    std::cerr << (*j).first.first << "-" << (*j).second.first << " "
         << (*j).first.second << "-" << (*j).second.second << std::endl;
    j++;
  }
  // Interval b(Key(320, 900),Key(330.1,910.7));
  Interval b(Key(320, 900),Key(330,910.7));
  Segment_tree_2.enclosing_query(b,std::back_inserter(N));
  j = N.begin();
  std::cerr << "\n enclosing_query (320, 900),(330.1,910.7)\n";
  while(j!=N.end())
  {
    std::cerr << (*j).first.first << "-" << (*j).second.first << " "
         << (*j).first.second << "-" << (*j).second.second << std::endl;
    j++;
  }
  if(Segment_tree_2.segment_tree_2->is_valid())
    std::cerr << "Tree  is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0;
}
