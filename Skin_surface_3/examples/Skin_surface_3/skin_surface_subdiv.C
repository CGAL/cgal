// examples/Skin_surface_3/skin_surface_subdiv.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <list>

#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Skin_surface_3::Bare_point                          Bare_point;

// Each facet has a pointer to the tetrahedron of the TMC it is contained in
typedef Skin_surface_3::Triangulated_mixed_complex          TMC;
typedef CGAL::Skin_surface_polyhedral_items_3<TMC>          Poly_items;
typedef CGAL::Polyhedron_3<K,Poly_items>                    Polyhedron;

int main(int argc, char *argv[]) {
  std::list<Weighted_point> l;
  FT                        shrinkfactor = 0.5;

  l.push_front(Weighted_point(Bare_point(0,0,0), 1));
  l.push_front(Weighted_point(Bare_point(0,1,0), 2));
  l.push_front(Weighted_point(Bare_point(0,0,2), 1));

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  CGAL::subdivide_skin_surface_mesh_3(p, skin_surface);

  std::ofstream out("mesh.off");
  out << p;

  return 0;
}
