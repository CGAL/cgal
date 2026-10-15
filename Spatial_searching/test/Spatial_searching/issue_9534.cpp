#include <CGAL/Dimension.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_d.h>

#include <array>
#include <iterator>
#include <vector>

template<class Kernel>
int test()
{
  using Traits = CGAL::Search_traits_d<Kernel>;
  using Point = typename Traits::Point_d;
  using FT = typename Kernel::FT;
  using Tree = CGAL::Kd_tree<Traits>;
  using Fuzzy_iso_box = CGAL::Fuzzy_iso_box<Traits>;
  using Fuzzy_sphere = CGAL::Fuzzy_sphere<Traits>;

  const std::array<FT, 2> c0 = {FT(0), FT(0)};
  const std::array<FT, 2> c1 = {FT(1), FT(1)};
  const std::array<FT, 2> c2 = {FT(2), FT(2)};

  const Point p0(2, c0.begin(), c0.end());
  const Point p1(2, c1.begin(), c1.end());
  const Point p2(2, c2.begin(), c2.end());

  Tree tree;
  tree.insert(p0);
  tree.insert(p1);
  tree.insert(p2);

  const Fuzzy_iso_box box(p0, p2);
  std::vector<Point> result;
  tree.search(std::back_inserter(result), box);

  if (result.empty())
    return 1;

  result.clear();

  const Fuzzy_sphere sphere(p1, FT(2));
  tree.search(std::back_inserter(result), sphere);

  if (result.empty())
    return 2;

  tree.remove(p1);

  if (tree.size() != 2)
    return 3;

  return 0;
}

int main()
{
  using Static_kernel = CGAL::Epeck_d<CGAL::Dimension_tag<2> >;
  using Dynamic_kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;

  int result = test<Static_kernel>();

  if (result != 0)
    return result;

  result = test<Dynamic_kernel>();

  if (result != 0)
    return 10 + result;

  return 0;
}
