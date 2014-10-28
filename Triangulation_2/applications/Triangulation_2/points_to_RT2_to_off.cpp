#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_filtered_traits_2<K>  Traits;
typedef CGAL::Regular_triangulation_2<Traits> Regular_triangulation;

int main()
{
  std::ifstream in("data/points.cin");

  Regular_triangulation::Weighted_point wp;
  std::vector<Regular_triangulation::Weighted_point> wpoints;

  while(in >> wp)
    wpoints.push_back(wp);

  Regular_triangulation rt(wpoints.begin(), wpoints.end());
  CGAL_assertion(rt.is_valid(true));
  std::ofstream off_stream("data/rt2.off");
  CGAL::export_triangulation_2_to_off(off_stream, rt);
  return 0;
}
