#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <variant>

template <class K>
void check()
{
  using Point = typename K::Point_3;
  using Triangle = typename K::Triangle_3;
  using Segment = typename K::Segment_3;
  const Point a(0, 1.0208333333314730e-1, 2.9199999999999998e-1);
  const Point b(0, 8.7499999999775646e-2, 2.9199999999999998e-1);
  const Point c(1.0617182181505985e-2, 8.7299626865148805e-2, 2.9199898847343825e-1);
  const Point d(0, 7.2916666666482610e-2, 2.9199999999999998e-1);
  const std::array<Point, 3> first = {a, b, c}, second = {d, b, c};
  std::array<int, 3> i = {0, 1, 2};
  do
  {
    std::array<int, 3> j = {0, 1, 2};
    do
    {
      const Triangle t1(first[i[0]], first[i[1]], first[i[2]]);
      const Triangle t2(second[j[0]], second[j[1]], second[j[2]]);
      for(bool reverse : {false, true})
      {
        const auto result = reverse ? CGAL::intersection(t2, t1) : CGAL::intersection(t1, t2);
        assert(result);
        const Segment* segment = std::get_if<Segment>(&*result);
        assert(segment);
        assert((segment->source() == b && segment->target() == c) ||
               (segment->source() == c && segment->target() == b));
      }
    } while(std::next_permutation(j.begin(), j.end()));
  } while(std::next_permutation(i.begin(), i.end()));
}

int main()
{
  check<CGAL::Exact_predicates_inexact_constructions_kernel>();
  check<CGAL::Exact_predicates_exact_constructions_kernel>();
}
