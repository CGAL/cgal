#include <CGAL/Interval_skip_list.h>
#include <CGAL/Interval_skip_list_interval.h>
#include <vector>
#include <iostream>
#include <list>
typedef CGAL::Interval_skip_list_interval<double> Interval;
typedef CGAL::Interval_skip_list<Interval> Interval_skip_list;

void
fct()
{
 
  Interval_skip_list isl;
  int i, n, d;

  n = 10;
  d = 3;
  //std::cin >> n >> d;
  std::vector<Interval> intervals(n);
  for(i = 0; i < n; i++) {
    intervals[i] = Interval(i,i+d);
  }
  std::random_shuffle(intervals.begin(), intervals.end());

  for(i = 0; i < n; i++) {
    isl.insert(intervals[i]);
  }

  for(i = 0; i < n+d; i++) {
    std::list<Interval> L;
    isl.find_intervals(i, std::back_inserter(L));
    for(std::list<Interval>::iterator it = L.begin(); it != L.end(); it++){
      std::cout << *it;
    }  
    std::cout << std::endl;
  }

  std::cout << isl;
  std::cout << std::endl;

  std::cout << * isl.begin();
  for(i = 0; i < n; i++) {
    isl.remove(intervals[i]);
  }
  std::cout << isl;
  std::cout << std::endl;
  
}

int
main()
{

  fct();
  return 0;

}

  
