// examples/Skin_surface_3/skin_surface_simple.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>

#include <list>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::RT                                  RT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;

bool test(char * filename, double shrink) {
  std::list<Weighted_point> l;
  std::ifstream in(filename);
  CGAL_assertion(in.is_open());
  Weighted_point wp;
  while (in >> wp) l.push_front(wp);

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrink);

  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  return (p.is_valid() && p.is_closed());
}

int main(int argc, char *argv[]) {
  bool result;
  char *filename;

  filename = "data/caffeine.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/ball.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/degenerate.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test1.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test2.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test3.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test4.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test5.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test6.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test7.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test8.cin";
  CGAL_assertion(result);
  result = test(filename, .5);
  filename = "data/test9.cin";
  CGAL_assertion(result);
  result = test(filename, .5);
  filename = "data/test10.cin";
  result = test(filename, .5);
  CGAL_assertion(result);
  filename = "data/test11.cin";
  result = test(filename, .5);
  CGAL_assertion(result);

  filename = "data/caffeine.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/ball.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/degenerate.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test1.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test2.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test3.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test4.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test5.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test6.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test7.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test8.cin";
  CGAL_assertion(result);
  result = test(filename, .85);
  filename = "data/test9.cin";
  CGAL_assertion(result);
  result = test(filename, .85);
  filename = "data/test10.cin";
  result = test(filename, .85);
  CGAL_assertion(result);
  filename = "data/test11.cin";
  result = test(filename, .85);
  CGAL_assertion(result);

  return 0;
}
