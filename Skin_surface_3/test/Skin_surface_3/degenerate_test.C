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
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Skin_surface_3::Bare_point                          Bare_point;
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

  result = test("data/caffeine.cin", .5);
  CGAL_assertion(result);
  result = test("data/ball.cin", .5);
  CGAL_assertion(result);
  result = test("data/degenerate.cin", .5);
  CGAL_assertion(result);
  result = test("data/test1.cin", .5);
  CGAL_assertion(result);
  result = test("data/test2.cin", .5);
  CGAL_assertion(result);
  result = test("data/test3.cin", .5);
  CGAL_assertion(result);
  result = test("data/test4.cin", .5);
  CGAL_assertion(result);
  result = test("data/test5.cin", .5);
  CGAL_assertion(result);
  result = test("data/test6.cin", .5);
  CGAL_assertion(result);
  result = test("data/test7.cin", .5);
  CGAL_assertion(result);
  result = test("data/test8.cin", .5);
  CGAL_assertion(result);
  result = test("data/test9.cin", .5);
  CGAL_assertion(result);
  result = test("data/test10.cin", .5);
  CGAL_assertion(result);
  result = test("data/test11.cin", .5);
  CGAL_assertion(result);

  result = test("data/caffeine.cin", .85);
  CGAL_assertion(result);
  result = test("data/ball.cin", .85);
  CGAL_assertion(result);
  result = test("data/degenerate.cin", .85);
  CGAL_assertion(result);
  result = test("data/test1.cin", .85);
  CGAL_assertion(result);
  result = test("data/test2.cin", .85);
  CGAL_assertion(result);
  result = test("data/test3.cin", .85);
  CGAL_assertion(result);
  result = test("data/test4.cin", .85);
  CGAL_assertion(result);
  result = test("data/test5.cin", .85);
  CGAL_assertion(result);
  result = test("data/test6.cin", .85);
  CGAL_assertion(result);
  result = test("data/test7.cin", .85);
  CGAL_assertion(result);
  result = test("data/test8.cin", .85);
  CGAL_assertion(result);
  result = test("data/test9.cin", .85);
  CGAL_assertion(result);
  result = test("data/test10.cin", .85);
  CGAL_assertion(result);
  result = test("data/test11.cin", .85);
  CGAL_assertion(result);

  return 0;
}
