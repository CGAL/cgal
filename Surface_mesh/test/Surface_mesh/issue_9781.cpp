#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <cassert>
#include <utility>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;

int main()
{
  Mesh source;
  const Mesh::Vertex_index v0 = source.add_vertex(Point(0, 1, 2));
  const Mesh::Vertex_index v1 = source.add_vertex(Point(3, 4, 5));

  Mesh destination;
  destination = std::move(source);
  assert(source.is_valid());
  assert(destination.is_valid());
  assert(source.is_empty());
  assert((destination.point(v0) == Point(0, 1, 2)));
  assert((destination.point(v1) == Point(3, 4, 5)));

  const Mesh::Vertex_index reused = source.add_vertex(Point(6, 7, 8));
  assert((destination.point(v0) == Point(0, 1, 2)));
  assert((destination.point(v1) == Point(3, 4, 5)));

  Mesh constructed(std::move(source));
  assert(source.is_valid());
  assert(constructed.is_valid());
  assert(source.is_empty());
  assert((constructed.point(reused) == Point(6, 7, 8)));

  return 0;
}
