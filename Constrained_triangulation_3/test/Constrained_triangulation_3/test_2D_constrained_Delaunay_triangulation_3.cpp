#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using PolygonMesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = PolygonMesh::Vertex_index;

int main()
{
  const std::array<Point, 4> points{ {
    {0., 0., 0.},
    {1., 0., 0.},
    {1., 1., 0.},
    {0., 1., 0.}
  } };
  const std::array<std::array<std::size_t, 4>, 1> polygons{ {
    {0, 1, 2, 3}
  } };
  std::array<vertex_descriptor, 4> vertices;
  PolygonMesh mesh;
  std::transform(points.begin(), points.end(), vertices.begin(),
    [&mesh](const Point& p) { return mesh.add_vertex(p); });

  [[maybe_unused]] auto fd = CGAL::Euler::add_face(vertices, mesh);
  assert(CGAL::is_valid_polygon_mesh(mesh));

  auto cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh);
  assert(cdt.is_valid());
  assert(cdt.constrained_facets().size() == 2);

  cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
  assert(cdt.is_valid());
  assert(cdt.constrained_facets().size() == 2);
}
