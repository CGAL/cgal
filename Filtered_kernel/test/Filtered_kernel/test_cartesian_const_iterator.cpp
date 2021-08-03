#include <algorithm>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using Traits_2 = CGAL::Search_traits_2<Kernel>;
using Traits_3 = CGAL::Search_traits_3<Kernel>;

int main() {

  // Testing 2D.
  Traits_2 traits_2;
  auto construct_it_2 = traits_2.construct_cartesian_const_iterator_d_object();
  const Point_2 pp2(1,2);
  const Point_2 qq2(2,3);
  const auto p2 = construct_it_2(pp2);
  const auto q2 = construct_it_2(qq2);
  const FT d2 = p2[0] - q2[0];
  assert(d2 < FT(0));

  const auto plen2 = std::distance(construct_it_2(pp2), construct_it_2(pp2, 0));
  const auto qlen2 = std::distance(construct_it_2(qq2), construct_it_2(qq2, 0));
  assert(plen2 == 2 && qlen2 == 2);

  // Testing 3D.
  Traits_3 traits_3;
  auto construct_it_3 = traits_3.construct_cartesian_const_iterator_d_object();
  const Point_3 pp3(1,2,3);
  const Point_3 qq3(2,3,4);
  const auto p3 = construct_it_3(pp3);
  const auto q3 = construct_it_3(qq3);
  const FT d3 = q3[2] - p3[2];
  assert(d3 > FT(0));

  const auto plen3 = std::distance(construct_it_3(pp3), construct_it_3(pp3, 0));
  const auto qlen3 = std::distance(construct_it_3(qq3), construct_it_3(qq3, 0));
  assert(plen3 == 3 && qlen3 == 3);

  return EXIT_SUCCESS;
}
