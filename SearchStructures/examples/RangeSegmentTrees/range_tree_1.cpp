// Implementation: Testprogram for 1-dimensional Range Trees
// A two dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.

#include <CGAL/Range_tree_k.h>
#include "include/Tree_Traits.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <utility>



typedef  CGAL::Range_tree_1< CGAL::Tree_traits_1> Range_tree_1_type;

int main()
{
  typedef  CGAL::Tree_traits_1::Key Key;
  typedef  CGAL::Tree_traits_1::Interval Interval;
  std::vector<Key> InputList, OutputList;
  std::vector<Key>::iterator first, last, current;

  InputList.push_back(8.0);
  InputList.push_back(1.0);
  InputList.push_back(3.9);
  InputList.push_back(2.0);
  InputList.push_back(5.0);
  InputList.push_back(4.8);
  InputList.push_back(7.7);
  InputList.push_back(6.8);

  first = InputList.begin();
  last = InputList.end();

  Range_tree_1_type Range_tree_1(first,last);

  Interval win(Interval(4.6, 6.8));

  std::cerr << "\n Window Query (4.6, 6.8) \n";
  Range_tree_1.window_query(win, std::back_inserter(OutputList));
  current=OutputList.begin();

  while(current!=OutputList.end())
  {
    std::cerr << (*current) << std::endl;
    current++;
  }

  if(Range_tree_1.range_tree_1->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0;
}
