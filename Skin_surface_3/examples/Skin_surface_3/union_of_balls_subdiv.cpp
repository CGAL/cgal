#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Union_of_balls_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_union_of_balls_3.h>
#include <CGAL/subdivide_union_of_balls_mesh_3.h>
#include <list>

#include "skin_surface_writer.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Union_of_balls_3<Traits>                      Union_of_balls_3;
typedef Union_of_balls_3::Weighted_point                    Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K,
  CGAL::Skin_surface_polyhedral_items_3<Union_of_balls_3> >   Polyhedron;

int main() {
  std::list<Weighted_point> l;

  l.push_front(Weighted_point(Bare_point( 1,-1,-1), 1.25));
  l.push_front(Weighted_point(Bare_point( 1, 1, 1), 1.25));
  l.push_front(Weighted_point(Bare_point(-1, 1,-1), 1.25));
  l.push_front(Weighted_point(Bare_point(-1,-1, 1), 1.25));

  Polyhedron p;

  Union_of_balls_3 union_of_balls(l.begin(), l.end(), Traits(), true);
  CGAL::mesh_union_of_balls_3(union_of_balls, p);

  CGAL::subdivide_union_of_balls_mesh_3(union_of_balls, p);

  std::ofstream out("output.off");
  out << p;
  return 0;
}
