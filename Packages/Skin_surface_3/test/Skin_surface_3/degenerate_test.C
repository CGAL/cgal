// test/Skin_surface_3/degenerate_test.C
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>

#include <list>
#include <fstream>

typedef CGAL::Skin_surface_traits_3<>                    Skin_surface_traits;
typedef Skin_surface_traits::Regular_traits              Regular_traits;
typedef Regular_traits::Bare_point                       Reg_point;
typedef Regular_traits::Weighted_point                   Reg_weighted_point;
typedef CGAL::Polyhedron_3<Skin_surface_traits::Polyhedron_kernel> Polyhedron;

bool test(char * filename, double shrink) {
  std::list<Reg_weighted_point> l;
  std::ifstream in(filename);
  CGAL_assertion(in.is_open());
  Reg_weighted_point wp;
  while (in >> wp) l.push_front(wp);

  Polyhedron p;
  Skin_surface_traits skin_surface_traits(shrink);
  CGAL::skin_surface_3(l.begin(), l.end(), p, skin_surface_traits);
  
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

  result = test("data/caffeine.cin", .85);
  CGAL_assertion(result);
  result = test("data/ball.cin", .85);
  CGAL_assertion(result);
  result = test("data/degenerate.cin", .85);
  CGAL_assertion(result);

  return 0;
}
