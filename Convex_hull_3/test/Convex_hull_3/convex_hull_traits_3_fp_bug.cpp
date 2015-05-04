#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <fstream>
#include <cassert>

typedef CGAL::Epick K;

int main()
{
  K::Point_3 p;
  std::set<K::Point_3> pointset;
  std::ifstream input("convex_hull_traits_3_fp_bug.xyz");

  assert(input);
  while ( input >> p )
    pointset.insert(p);

  CGAL::Polyhedron_3<K> r;
  CGAL::convex_hull_3(pointset.begin(), pointset.end(), r);

  assert(r.size_of_vertices()==82);

  return 0;
}
