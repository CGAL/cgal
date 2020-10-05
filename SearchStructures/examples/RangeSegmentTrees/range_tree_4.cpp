// Implementation: Testprogram for 4-dimensional Range Trees
// A four dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.

#include <CGAL/Range_tree_k.h>
#include "include/Tree_Traits.h"
#include <iostream>
#include <list>
#include <iterator>

typedef CGAL::Range_tree_4<CGAL::Tree_traits_4> Range_tree_4_type;

int main()
{
  typedef CGAL::Tree_traits_4::Key Key;
  typedef CGAL::Tree_traits_4::Interval Interval;

  std::list<Key> InputList, OutputList;
  std::list<Key>::iterator first, last, current;

  InputList.push_back(Key(8,5.1,34, 1.11));
  InputList.push_back(Key(1,1.1,67, 1.23));
  InputList.push_back(Key(3,2.1,56, 1.34));
  InputList.push_back(Key(2,6.1,89, 1.09));
  InputList.push_back(Key(5,4.1,45, 1.009));
  InputList.push_back(Key(4,8.1,76, 1.98));
  InputList.push_back(Key(7,7.1,87, 1.333));
  InputList.push_back(Key(6,3.1,78, 1.45));

  first = InputList.begin();
  last = InputList.end();

  Range_tree_4_type Range_tree_4(first,last);

  Key a=Key(4,8.1,45, 1.12);
  Key b=Key(5,8.2,89, 1.99);
  Interval win=Interval(a,b);

  std::cerr << "\n Window Query:(4,8.1,45, 1.12), (5,8.2,89, 1.99)\n";
  Range_tree_4.window_query(win, std::back_inserter(OutputList));
  current=OutputList.begin();

  while(current!=OutputList.end())
  {
    std::cerr << (*current).key_1<< "-" <<  (*current).key_2 << "-"
         <<  (*current).key_3 << (*current).key_4 << std::endl;
    current++;
  }

  if(Range_tree_4.range_tree_4->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";

  return 0;
}
