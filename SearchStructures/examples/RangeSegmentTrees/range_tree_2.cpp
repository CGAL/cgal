// Implementation: Testprogram for 2-dimensional Range Trees
// A two dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.

#include <CGAL/Range_tree_k.h>
#include "include/Tree_Traits.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <utility>

typedef CGAL::Range_tree_2<CGAL::Tree_traits_2> Range_tree_2_type;

int main()
{
  typedef CGAL::Tree_traits_2::Key Key;
  typedef CGAL::Tree_traits_2::Interval Interval;

  std::vector<Key> InputList, OutputList;
  typedef std::vector<Key>::iterator V_iterator;
  V_iterator first, last, current;

  int i,j;
  j= 100;
  int first_key=1, second_key=3;
  for(i=1;i<j;i++)
  {
    InputList.push_back(Key(first_key++,second_key++));
  }

  first = InputList.begin();
  last = InputList.end();
  std::cerr << InputList.size();

  Range_tree_2_type *Range_tree_2 = new Range_tree_2_type(first,last);

  Key a=Key(4,8.1);
  Key b=Key(5,8.2);
  Interval win=Interval(a,b);

  std::cerr << "\n Window Query: \n";
  Range_tree_2->window_query(win, std::back_inserter(OutputList));
  current=OutputList.begin();

  while(current!=OutputList.end())
  {
    std::cerr << (*current).first<< "-" << (*current).second << std::endl;
    current++;
  }

  if(Range_tree_2->range_tree_2->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";

  return 0;
}
