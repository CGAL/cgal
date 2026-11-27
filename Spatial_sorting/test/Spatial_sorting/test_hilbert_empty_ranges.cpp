#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

using K = CGAL::Simple_cartesian<double>;
using Kd = CGAL::Cartesian_d<double>;
using Policy = CGAL::Hilbert_sort_middle_policy;

template <class Kernel, class Sort, class Point, class SortPolicy>
void check(const Kernel& kernel, const Sort& sort, const Point& point, SortPolicy policy)
{
  std::vector<Point> points;
  sort(points.begin(), points.end());
  CGAL::hilbert_sort(points.begin(), points.end(), kernel, policy);
  CGAL::spatial_sort(points.begin(), points.end(), kernel, policy);
  CGAL::spatial_sort(points.begin(), points.end(), kernel);
  assert(points.empty());
  for(std::size_t size : {1, 2, 17}) {
    points.assign(size, point);
    const auto original = points;
    sort(points.end(), points.end());
    CGAL::hilbert_sort(points.end(), points.end(), kernel, policy);
    CGAL::spatial_sort(points.end(), points.end(), kernel, policy);
    CGAL::spatial_sort(points.end(), points.end(), kernel);
    assert(points == original);
    sort(points.begin(), points.end());
    CGAL::hilbert_sort(points.begin(), points.end(), kernel, policy);
    CGAL::spatial_sort(points.begin(), points.end(), kernel, policy);
    CGAL::spatial_sort(points.begin(), points.end(), kernel);
    assert(points == original);
  }
}

int main(int argc, char** argv)
{
  const std::string dimension = argc > 1 ? argv[1] : "all";
  if(dimension == "all" || dimension == "2")
    check(K(), CGAL::Hilbert_sort_2<K, Policy>(), K::Point_2(0, 1), Policy());
  if(dimension == "all" || dimension == "3")
    check(K(), CGAL::Hilbert_sort_3<K, Policy>(), K::Point_3(0, 1, 2), Policy());
  if(dimension == "all" || dimension == "d") {
    const double coordinates[] = {0, 1, 2, 3};
    check(Kd(), CGAL::Hilbert_sort_d<Kd, Policy>(),
          Kd::Point_d(4, coordinates, coordinates + 4), Policy());
  }
  if(dimension == "all" || dimension == "median-d") {
    const double coordinates[] = {0, 1, 2, 3};
    using Median = CGAL::Hilbert_sort_median_policy;
    check(Kd(), CGAL::Hilbert_sort_d<Kd, Median>(),
          Kd::Point_d(4, coordinates, coordinates + 4), Median());
  }
  std::cout << "Empty-range regression passed (" << dimension << ").\n";
}
