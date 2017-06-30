#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/IO/Triangulation_off_ostream_3.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef CGAL::Regular_triangulation_3<Traits> Regular_triangulation;

int main()
{
  std::ifstream in("data/points.cin");

  Regular_triangulation::Weighted_point wp;
  std::vector<Regular_triangulation::Weighted_point> wpoints;

  while(in >> wp)
    wpoints.push_back(wp);

  Regular_triangulation rt(wpoints.begin(), wpoints.end());
  std::ofstream off_stream("data/rt3.off");
  CGAL::export_triangulation_3_to_off(off_stream, rt);
  return 0;
}
