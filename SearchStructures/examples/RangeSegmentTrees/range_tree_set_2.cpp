// Implementation: Testprogram for 2-dimensional Range Trees
// A two dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>
#include <iostream>
#include <utility>
#include <vector>
#include <iterator>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Range_segment_tree_set_traits_2<K> Traits;
typedef CGAL::Range_tree_2<Traits> Range_tree_2_type;

int main()
{
  typedef Traits::Key Key;
  typedef Traits::Interval Interval;
  //  Range_tree_2_type RR;
  std::vector<Key> InputList, OutputList;
  std::vector<Key>::iterator first, last, current;

  InputList.push_back(Key(8,5.1));
  InputList.push_back(Key(1,1.1));
  InputList.push_back(Key(3,2.1));
  InputList.push_back(Key(2,6.1));
  InputList.push_back(Key(5,4.1));
  InputList.push_back(Key(4,8.1));
  InputList.push_back(Key(7,7.1));
  InputList.push_back(Key(6,3.1));

  first = InputList.begin();
  last = InputList.end();

  Range_tree_2_type Range_tree_2(first,last);

  Key a=Key(4,8.1);
  Key b=Key(5,8.2);
  Interval win=Interval(a,b);


  std::cerr << std::endl << "Window Query: lower left point: (4.0,5.0),";
  std::cerr << "upper right point: (8.1,8.2)" << std::endl;
  Range_tree_2.window_query(win, std::back_inserter(OutputList));
  current=OutputList.begin();

  while(current!=OutputList.end())
  {
    std::cerr << (*current).x()<< "-" << (*current).y() << std::endl;
    current++;
  }

  if(Range_tree_2.range_tree_2->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";

  return 0;
}
