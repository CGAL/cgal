#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <list>

#include <fstream>
#include "skin_surface_writer.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Skin_surface_traits_3<K>                        Traits;
typedef CGAL::Skin_surface_3<Traits>                          Skin_surface_3;
typedef Skin_surface_3::FT                                    FT;
typedef Skin_surface_3::Weighted_point                        Weighted_point;
typedef Weighted_point::Point                                 Bare_point;
typedef CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> Polyhedral_items;
typedef CGAL::Polyhedron_3<K, Polyhedral_items>               Polyhedron;

int main() {
  std::list<Weighted_point> l;
  FT                        shrinkfactor = 0.5;

  l.push_front(Weighted_point(Bare_point(0,0,0), 1));
  l.push_front(Weighted_point(Bare_point(0,1,0), 2));
  l.push_front(Weighted_point(Bare_point(0,0,2), 1));

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  CGAL::subdivide_skin_surface_mesh_3(skin_surface, p);

  std::ofstream out("mesh.off");
  write_polyhedron_with_normals(skin_surface, p, out);

  return 0;
}
