#include <CGAL/algorithm.h>
#include <vector>
#include <iostream>

using std::vector;
using std::pair;
using std::cout;
using std::endl;
using CGAL::min_max_element;

int main()
{
  vector< int > v;
  v.push_back(3);
  v.push_back(6);
  v.push_back(5);
  typedef std::vector< int >::iterator iterator;
  pair< iterator, iterator > p = min_max_element(v.begin(), v.end());
  cout << "min = " << *p.first << ", max = " << *p.second << endl;
  return 0;
}
