// examples/Skin_surface_3/skin_surface_pdb_reader.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>

// NGHK: for subdivision
#include <CGAL/subdivide_skin_surface_mesh_3.h>

#include <extract_balls_from_pdb.h>

#include <list>
#include "skin_surface_writer.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef CGAL::Simple_cartesian<float>                       Poly_K;
typedef Skin_surface_3::Triangulated_mixed_complex          TMC;
typedef CGAL::Skin_surface_polyhedral_items_3<TMC>          Poly_items;
typedef CGAL::Polyhedron_3<Poly_K,Poly_items>               Polyhedron;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <pdb-file>" << std::endl;
    return 1;
  }

  std::list<Weighted_point> l;
  FT                        shrinkfactor = 0.5;

  extract_balls_from_pdb(argv[1], K(), std::back_inserter(l));
  std::cout << "Read pdb" << std::endl;

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor, true, Traits(), true);
  std::cout << "Constructed Skin_surface_3" << std::endl;

  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);
  std::cout << "Meshed Skin_surface_3" << std::endl;

  //  CGAL::subdivide_skin_surface_mesh_3(p, skin_surface);
  //  std::cout << "Subdivided Skin_surface_3" << std::endl;

  std::ofstream out("mesh.off");
  write_polyhedron_with_normals(p, skin_surface, out);
  //out << p; 

  return 0;
}
