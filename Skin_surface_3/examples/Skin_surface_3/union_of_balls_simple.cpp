#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Union_of_balls_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_union_of_balls_3.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Union_of_balls_3<Traits>                      Union_of_balls_3;

typedef Union_of_balls_3::Bare_point                        Bare_point;
typedef Union_of_balls_3::Weighted_point                    Weighted_point;

typedef CGAL::Polyhedron_3<K>                               Polyhedron;

int main()
{
  std::list<Weighted_point> l;

  l.push_front(Weighted_point(Bare_point(0,0,0), 1));
  l.push_front(Weighted_point(Bare_point(0,1,0), 2));
  l.push_front(Weighted_point(Bare_point(0,0,2), 1));

  Union_of_balls_3 union_of_balls(l.begin(), l.end());

  Polyhedron p;
  CGAL::mesh_union_of_balls_3(union_of_balls, p);

  return 0;
}
