#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Kd_tree.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using Traits = CGAL::Search_traits_2<Kernel>;
using Point = Kernel::Point_2;

template <typename Cache>
void check_search()
{
  using Tree = CGAL::Kd_tree<Traits, CGAL::Sliding_midpoint<Traits>,
                             CGAL::Tag_true, Cache>;
  const std::vector<Point> points = {
    Point(0, 0), Point(3, 4), Point(5, 0), Point(-5, 0),
    Point(6, 0), Point(4, 4)
  };
  const Tree tree(points.begin(), points.end());
  std::vector<Point> found;
  tree.search(std::back_inserter(found), CGAL::Fuzzy_sphere<Traits>(Point(0, 0), 5));
  assert(found.size() == 4);
  for(std::size_t i = 0; i < 4; ++i)
    assert(std::find(found.begin(), found.end(), points[i]) != found.end());

  found.clear();
  tree.search(std::back_inserter(found), CGAL::Fuzzy_sphere<Traits>(Point(0, 0), 0));
  assert(found.size() == 1);
  assert(found.front() == Point(0, 0));
}

int main()
{
  const Traits traits;
  const auto coordinates = traits.construct_cartesian_const_iterator_d_object();
  const CGAL::Fuzzy_sphere<Traits> sphere(Point(0, 0), 5);
  for(const Point& point : {Point(0, 0), Point(3, 4), Point(5, 0), Point(6, 0)})
    assert(sphere.contains_point_given_as_coordinates(coordinates(point), coordinates(point, 0))
           == sphere.contains(point));

  check_search<CGAL::Tag_false>();
  check_search<CGAL::Tag_true>();
}
