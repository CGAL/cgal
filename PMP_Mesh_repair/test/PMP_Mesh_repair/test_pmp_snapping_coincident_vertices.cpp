#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/property_map.h>

#include <cassert>
#include <cstddef>
#include <initializer_list>

namespace PMP = CGAL::Polygon_mesh_processing;

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template <typename Mesh>
void add_triangle(Mesh& mesh, typename Mesh::Vertex_index a,
                 typename Mesh::Vertex_index b, typename Mesh::Vertex_index c,
                 bool reverse)
{
  auto face = reverse ? mesh.add_face(a, c, b) : mesh.add_face(a, b, c);
  assert(face != Mesh::null_face());
}

template <typename K>
void test_coincident_vertices(bool self_snapping, bool reverse, double height,
                             double tolerance)
{
  using Point = typename K::Point_3;
  using FT = typename K::FT;
  using Mesh = CGAL::Surface_mesh<Point>;

  Mesh source, other;
  Mesh& target = self_snapping ? source : other;
  const Point position(0, 0, height);
  auto p = source.add_vertex(position);
  auto a = source.add_vertex(Point(-1, 0, 0));
  auto b = source.add_vertex(Point(0, 0, 0));
  add_triangle(source, p, a, b, reverse);
  auto q = source.add_vertex(position);
  auto c = source.add_vertex(Point(0, 0, 30));
  auto d = source.add_vertex(Point(-1, 0, 0));
  add_triangle(source, q, c, d, reverse);
  auto e = target.add_vertex(Point(0, 0, 0));
  auto f = target.add_vertex(Point(0, 0, 30));
  auto g = target.add_vertex(Point(1, 0, 0));
  add_triangle(target, e, g, f, reverse);

  CGAL::Constant_property_map<typename Mesh::Vertex_index, FT> tol{FT(tolerance)};
  if(self_snapping)
    PMP::experimental::snap_borders(source, tol);
  else
    PMP::experimental::snap_borders(source, tol, target, tol);

  assert(source.point(p) == position);
  assert(source.point(q) == position);
  assert(CGAL::is_valid_polygon_mesh(source));
  assert(CGAL::is_valid_polygon_mesh(target));
  if(!self_snapping)
    CGAL::copy_face_graph(target, source);

  assert(source.number_of_vertices() == 10);
  const auto stitched = PMP::stitch_borders(source);
  assert(stitched == 3);
  assert(source.number_of_vertices() == 5);
  assert(source.number_of_faces() == 4);
  assert(CGAL::is_valid_polygon_mesh(source));

  std::size_t border_edges = 0;
  for(auto h : halfedges(source))
    if(is_border(h, source))
      ++border_edges;
  assert(border_edges == 4);
  auto components = source.template add_property_map<typename Mesh::Face_index,
                                                      std::size_t>("f:component").first;
  assert(PMP::connected_components(source, components) == 1);
}

template <typename K>
void test_projection()
{
  using Point = typename K::Point_3;
  using FT = typename K::FT;
  using Mesh = CGAL::Surface_mesh<Point>;

  Mesh source, target;
  const Point position(0.00001, 0, 0.1);
  auto p = source.add_vertex(position);
  auto a = source.add_vertex(Point(-1, 0, 0));
  auto b = source.add_vertex(Point(-1, 0, 1));
  add_triangle(source, p, a, b, false);
  auto e = target.add_vertex(Point(0, 0, 0));
  auto f = target.add_vertex(Point(0, 0, 30));
  auto g = target.add_vertex(Point(1, 1, 0));
  add_triangle(target, e, f, g, false);

  CGAL::Constant_property_map<typename Mesh::Vertex_index, FT> tol{FT(0.001)};
  PMP::experimental::snap_borders(source, tol, target, tol,
                                 CGAL::parameters::default_values(),
                                 CGAL::parameters::do_lock_mesh(true));
  assert(source.point(p) != position);
  assert(source.point(p).x() == 0);
  assert(source.point(p).y() == 0);
  assert(K().has_on_3_object()(typename K::Segment_3(Point(0, 0, 0), Point(0, 0, 30)),
                             source.point(p)));
  assert(target.point(e) == Point(0, 0, 0));
  assert(target.point(f) == Point(0, 0, 30));
  assert(target.point(g) == Point(1, 1, 0));
  assert(target.number_of_vertices() == 4);
  assert(CGAL::is_valid_polygon_mesh(source));
  assert(CGAL::is_valid_polygon_mesh(target));
}

template <typename K>
void test()
{
  for(bool self_snapping : {true, false})
    for(bool reverse : {false, true})
    {
      test_coincident_vertices<K>(self_snapping, reverse, 0.1, 1e-6);
      test_coincident_vertices<K>(self_snapping, reverse, 1.7347235e-14, 1e-16);
    }
  test_projection<K>();
}

int main()
{
  test<EPICK>();
  test<EPECK>();
}
