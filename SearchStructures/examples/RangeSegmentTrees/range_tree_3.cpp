// Implementation: Testprogram for 3-dimensional Range Trees
// A three dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.

#include <CGAL/Range_tree_k.h>
#include "include/Tree_Traits.h"
#include <iostream>
#include <list>
#include <iterator>
#include <utility>

typedef CGAL::Range_tree_3<CGAL::Tree_traits_3> Range_tree_3_type;


int main()
{
  typedef CGAL::Tree_traits_3::Key Key;
  typedef CGAL::Tree_traits_3::Interval Interval;

  std::list<Key> InputList, OutputList;
  std::list<Key>::iterator first, last, current;

  InputList.push_back(Key(8,5.1,34));
  InputList.push_back(Key(1,1.1,67));
  InputList.push_back(Key(3,2.1,56));
  InputList.push_back(Key(2,6.1,89));
  InputList.push_back(Key(5,4.1,45));
  InputList.push_back(Key(4,8.1,76));
  InputList.push_back(Key(7,7.1,87));
  InputList.push_back(Key(6,3.1,78));

  first = InputList.begin();
  last = InputList.end();

  Range_tree_3_type Range_tree_3(first,last);
  Key a=Key(4,8.1,45);
  Key b=Key(5,8.2,89);
  Interval win=Interval(a,b);

  std::cerr << "\n Window Query:(4,8.1,45),(5,8.2,89)\n";
  Range_tree_3.window_query(win, std::back_inserter(OutputList));
  current=OutputList.begin();

  while(current!=OutputList.end())
  {
    std::cerr << (*current).key_1<< "-" <<  (*current).key_2<< "-"
         << (*current).key_3 << std::endl;
    current++;
  }

  if(Range_tree_3.range_tree_3->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";

  return 0;
}
