// interval skip list test program  
// Author:  Eric N. Hanson, cis.ufl.edu, University of Florida

#include <string.h>

#include <CGAL/Interval_skip_list.h>
#include <vector>
#include <iostream>

#define maxKey 10000
#define offset 1000
#define maxIntervals 10000

void
fct()
{
  CGAL::IntervalSkipList isl;
  int i, n, d;

  std::cin >> n >> d;
  std::vector<CGAL::Interval> intervals(n);
  for(int i = 0; i < n; i++) {
    intervals[i] = CGAL::Interval('[', i, i+d, ']');
  }
  std::random_shuffle(intervals.begin(), intervals.end());

  for(i = 0; i < n; i++) {
    isl.insert(intervals[i]);
  }

  for(i = 0; i < n+d; i++) {
    std::list<CGAL::Interval> L;
    isl.findIntervals(i, std::back_inserter(L));
    for(std::list<CGAL::Interval>::iterator it = L.begin(); it != L.end(); it++){
      std::cout << *it;
    }  
    std::cout << std::endl;
  }

  isl.print();
  std::cout << std::endl;

  for(i = 0; i < n; i++) {
    isl.remove(intervals[i]);
  }
  isl.print();
  std::cout << std::endl;
  
}

int
main()
{

  fct();
  return 0;

}


