#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/make_skin_surface_mesh_3.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3                                          Bare_point;
typedef CGAL::Weighted_point<Bare_point,K::RT>              Weighted_point;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;

int main() {
  std::list<Weighted_point> l;
  double                    shrinkfactor = 0.5;

  l.push_front(Weighted_point(Bare_point( 1,-1,-1), 1.25));
  l.push_front(Weighted_point(Bare_point( 1, 1, 1), 1.25));
  l.push_front(Weighted_point(Bare_point(-1, 1,-1), 1.25));
  l.push_front(Weighted_point(Bare_point(-1,-1, 1), 1.25));

  Polyhedron p;

  CGAL::make_skin_surface_mesh_3(p, l.begin(), l.end(), shrinkfactor);

  return 0;
}
