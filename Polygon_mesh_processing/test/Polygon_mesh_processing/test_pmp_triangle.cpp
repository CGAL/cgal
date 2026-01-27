#include <CGAL/Polygon_mesh_processing/triangle.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/generators.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Triangle_3 Triangle_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

#include <iostream>

int main()
{
  Surface_mesh tri;
  Point_3 p(0, 0, 0), q(1, 0, 0), r(0, 1,0);
  CGAL::make_triangle(p, q, r, tri);
  boost::graph_traits<Surface_mesh>::face_descriptor fd = *faces(tri).begin();
  Triangle_3 t = CGAL::Polygon_mesh_processing::triangle(fd, tri, CGAL::parameters::geom_traits(K()));

  assert((t[0] == p && t[1] == q &&  t[2] == r) ||
         (t[1] == p && t[2] == q &&  t[0] == r) ||
         (t[2] == p && t[0] == q &&  t[1] == r));

  return 0;
}
