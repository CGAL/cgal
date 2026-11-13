#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/partition_is_valid_2.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <list>
#include <vector>

template <class K>
void check(const std::vector<std::array<int, 2>>& coordinates, std::size_t expected_parts)
{
  using Traits = CGAL::Partition_traits_2<K>;
  using Polygon = typename Traits::Polygon_2;
  std::vector<typename K::Point_2> points;
  for(const auto& p : coordinates)
    points.emplace_back(p[0], p[1]);
  for(std::size_t start = 0; start < points.size(); ++start)
  {
    std::list<Polygon> result;
    CGAL::optimal_convex_partition_2(points.begin(), points.end(), std::back_inserter(result));
    assert(CGAL::convex_partition_is_valid_2(points.begin(), points.end(), result.begin(), result.end()));
    assert(result.size() == expected_parts);
    std::rotate(points.begin(), points.begin() + 1, points.end());
  }
}

template <class K>
void test()
{
  check<K>({{-320, -448}, {-320, -384}, {-128, 24}, {-128, 256},
             {-320, 64}, {-320, -192}, {-512, -160}}, 2);
  check<K>({{-640, -160}, {-640, -384}, {-448, -384}, {-448, -448},
             {-320, -448}, {-320, -384}, {-128, -384}, {-128, -152},
             {-128, -128}, {-128, 0}, {-128, 24}, {-128, 256},
             {-320, 256}, {-448, 64}, {-320, 64}, {-320, -192},
             {-448, -192}, {-512, -192}, {-512, -160}}, 5);
  check<K>({{0, 0}, {1, 0}, {2, 0}, {2, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}}, 1);
}

int main()
{
  test<CGAL::Exact_predicates_inexact_constructions_kernel>();
  test<CGAL::Exact_predicates_exact_constructions_kernel>();
}
