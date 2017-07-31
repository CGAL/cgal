#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Gm_polyhedron;

int main()
{
  Point_3 points[] = { Point_3(1.0, 0.0, 0.0), Point_3(0.0, 1.0, 0.0), Point_3(0.0, 0.0, 1.0), Point_3(0.0, 0.0, 0.0) };
  Gm_polyhedron P1;
  CGAL::convex_hull_3(points, &points[4], P1);
  CGAL_assertion( CGAL::Polygon_mesh_processing::is_outward_oriented(P1) );

  Point_3 points_bis[] = { Point_3(0.0, 1.0, 0.0), Point_3(1.0, 0.0, 0.0), Point_3(0.0, 0.0, 1.0), Point_3(0.0, 0.0, 0.0) };
  CGAL::convex_hull_3(points_bis, &points_bis[4], P1);
  CGAL_assertion( CGAL::Polygon_mesh_processing::is_outward_oriented(P1) );
}
