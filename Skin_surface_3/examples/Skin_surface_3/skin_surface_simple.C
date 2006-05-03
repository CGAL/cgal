// examples/Skin_surface_3/skin_surface_simple.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>

#include <list>
// NGHK: remove later
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Skin_surface_3::Bare_point                          Bare_point;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;

int main(int argc, char *argv[]) {
  std::list<Weighted_point> l;
  FT                        shrinkfactor = 0.5;

  Weighted_point wp;
  std::ifstream in("data/caffeine.cin");
  while (in >> wp) l.push_front(wp);
//   l.push_front(Weighted_point(Bare_point(0,0,0), 1));
//   l.push_front(Weighted_point(Bare_point(0,1,0), 2));
//   l.push_front(Weighted_point(Bare_point(0,0,2), 1));

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  std::ofstream out("mesh.off");
  out << p;

  return 0;
}
