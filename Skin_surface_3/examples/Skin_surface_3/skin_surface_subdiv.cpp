// examples/Skin_surface_3/skin_surface_subdiv.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mixed_complex_traits_3.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <list>


// NGHK: remove later
#include <fstream>
#include "skin_surface_writer.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mixed_complex_traits_3<K>                     Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::RT                                  RT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;

int main(int argc, char *argv[]) {
  std::list<Weighted_point> l;
  RT                        shrinkfactor = 0.5;

  l.push_front(Weighted_point(Bare_point(0,-1,-2), 1));
  l.push_front(Weighted_point(Bare_point(0,-1, 2), 1));
  l.push_front(Weighted_point(Bare_point(0, 1, 0), 1));

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);
  
  CGAL::subdivide_skin_surface_mesh_3(p, skin_surface);


  // NGHK: remove later
  std::ofstream out("mesh.off");
  write_polyhedron_with_normals(p, skin_surface, out);

  return 0;
}
