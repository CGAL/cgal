#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/convexity_check_3.h>
#include <fstream>
#include <cassert>

typedef CGAL::Epick K;
typedef CGAL::Epeck EK;

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

  CGAL::Polyhedron_3<EK> s;
  CGAL::copy_face_graph(r,s);
  assert(CGAL::is_strongly_convex_3(s));

  CGAL::Cartesian_converter<K, EK> to_EK;
  CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<EK>, EK> sotm(s);


  for(K::Point_3 p : pointset)
  {
    assert(sotm(to_EK(p)) != CGAL::ON_UNBOUNDED_SIDE);
  }

  return 0;
}
