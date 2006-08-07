// examples/Skin_surface_3/skin_surface_pdb_reader.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>

#include <extract_balls_from_pdb.h>

#include <list>
#include "skin_surface_writer.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Mixed_complex_traits_3<K>                       Traits;
typedef CGAL::Skin_surface_3<Traits>                          Skin_surface_3;
typedef Skin_surface_3::RT                                    RT;
typedef Skin_surface_3::Weighted_point                        Weighted_point;
typedef CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> Poly_items;
typedef CGAL::Polyhedron_3<K,Poly_items>                      Polyhedron;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <pdb-file>" << std::endl;
    return 0;
  }

  std::list<Weighted_point> l;
  RT                        shrinkfactor = 0.5;
  extract_balls_from_pdb(argv[1], K(), std::back_inserter(l));

  // Construct skin surface:
  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  // Extract mesh from the skin surface:
  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  // Subdivide skin surface
  // CGAL::subdivide_skin_surface_mesh_3(p, skin_surface, 1);

  std::ofstream out("mesh.off");
  write_polyhedron_with_normals(p, skin_surface, out);

  return 0;
}
