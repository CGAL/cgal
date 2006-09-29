// examples/Skin_surface_3/NGHK_skin_surface_subdiv.C
//#define CGAL_PROFILE
//#define CGAL_NO_ASSERTIONS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>
#include <list>

#include <fstream>
#include "skin_surface_writer.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::RT                                  RT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K,
  CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> >   Polyhedron;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <cin-file>" << std::endl;
    return 0;
  }

  std::list<Weighted_point> l;
  RT                        shrinkfactor = 0.5;

  Weighted_point wp;
  std::ifstream in(argv[1]);
  while (in >> wp) l.push_back(wp);

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor, false);

  Polyhedron p;

  std::cout << "Meshing ..." << std::endl;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  std::cout << "Subdividing ..." << std::endl;
  CGAL::subdivide_skin_surface_mesh_3(skin_surface, p, 1);

  std::cout << "Is closed: " << (p.is_closed() ? "Yes" : "No") << std::endl;

  std::ofstream out("mesh.off");
  write_polyhedron_with_normals(p, skin_surface, out);
//   out << p;

  return 0;
}
