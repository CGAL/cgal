#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_sphere_2.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>           Traits;
typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Traits> CDToS_2;

typedef CDToS_2::Point                                            Point;

int main()
{
  std::vector<K::Point_3> points;
  points.emplace_back( 2, 1, 1);
  points.emplace_back(-2, 1, 1);
  points.emplace_back( 1, 2, 1);
  points.emplace_back( 0, 1, 1);
  points.emplace_back( 1, 0, 1);
  points.emplace_back( 1, 1, 2);
  points.emplace_back( 1, 1, 0);

  Traits traits(Point(0, 0, 0));
  CDToS_2 dtos(traits);



}
