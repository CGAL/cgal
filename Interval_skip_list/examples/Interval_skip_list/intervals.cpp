#include <CGAL/Interval_skip_list.h>
#include <CGAL/Interval_skip_list_interval.h>
#include <CGAL/algorithm.h>
#include <vector>
#include <list>
#include <iostream>

typedef CGAL::Interval_skip_list_interval<double> Interval;
typedef CGAL::Interval_skip_list<Interval> Interval_skip_list;

int main()
{
  Interval_skip_list isl;
  int i, n, d;

  n = 10;
  d = 3;
  //std::cin >> n >> d;
  std::vector<Interval> intervals(n);
  for(i = 0; i < n; i++) {
    intervals[i] = Interval(i, i+d);
  }
  CGAL::cpp98::random_shuffle(intervals.begin(), intervals.end());

  isl.insert(intervals.begin(), intervals.end());

  for(i = 0; i < n+d; i++) {
    std::list<Interval> L;
    isl.find_intervals(i, std::back_inserter(L));
    for(std::list<Interval>::iterator it = L.begin(); it != L.end(); it++){
      std::cout << *it;
    }
    std::cout << std::endl;
  }

  for(i = 0; i < n; i++) {
    isl.remove(intervals[i]);
  }
  return 0;

}
