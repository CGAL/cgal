#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include "skin_surface_writer.h"
#include <list>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K,
  CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> >   Polyhedron;

int main() {
  std::list<Weighted_point> l;
  FT                        shrinkfactor = 0.5;

  l.push_front(Weighted_point(Bare_point( 1,-1,-1), 1.25));
  l.push_front(Weighted_point(Bare_point( 1, 1, 1), 1.25));
  l.push_front(Weighted_point(Bare_point(-1, 1,-1), 1.25));
  l.push_front(Weighted_point(Bare_point(-1,-1, 1), 1.25));

  Polyhedron p;

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);
  CGAL::mesh_skin_surface_3(skin_surface, p);

  CGAL::subdivide_skin_surface_mesh_3(skin_surface, p);

  std::ofstream out("mesh.off");
  out << p;

  return 0;
}
