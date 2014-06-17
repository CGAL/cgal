#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/IO/Triangulation_off_ostream.h>

#include <fstream>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef CGAL::Regular_triangulation_euclidean_traits<K>  Traits;
typedef CGAL::Regular_triangulation<Traits> RT;

int main()
{
  std::ifstream in("data/points.cin");

  RT::Weighted_point wp;
  std::vector<RT::Weighted_point> wpoints;

  int dim;
  in >> dim;
  while(in >> wp)
    wpoints.push_back(wp);

  // Build the Regular Triangulation
  RT rt(dim);
  rt.insert(wpoints.begin(), wpoints.end());
  std::ofstream off_stream("data/rt.off");
  CGAL::export_triangulation_to_off(off_stream, rt);
  return 0;
}
