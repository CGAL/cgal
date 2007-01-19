// examples/Skin_surface_3/skin_surface_pdb_reader.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K,
  CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> >   Polyhedron;

#include <extract_balls_from_pdb.h>

#include <list>
#include <fstream>
//#include "skin_surface_writer.h"


int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <pdb-file>" << std::endl;
    return 0;
  }

  std::list<Weighted_point> l;
  double                    shrinkfactor = 0.5;

  // Retrieve input balls:
  extract_balls_from_pdb(argv[1], K(), std::back_inserter(l));

  // Construct skin surface:
  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  // Extract mesh from the skin surface:
  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

//  std::ofstream out("mesh.off");
//  write_polyhedron_with_normals(p, skin_surface, out);

  return 0;
}
