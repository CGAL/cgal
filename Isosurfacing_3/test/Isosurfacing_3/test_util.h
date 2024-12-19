#ifndef CGAL_ISOSURFACING_TEST_UTIL_H
#define CGAL_ISOSURFACING_TEST_UTIL_H

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <cassert>
#include <iostream>

template <typename PointRange, typename PolygonRange>
bool has_duplicate_points(PointRange points, PolygonRange polygons) // intentional copies
{
  return CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, polygons) != 0;
}

template <typename PointRange, typename PolygonRange>
bool has_duplicate_polygons(PointRange points, PolygonRange polygons)
{
  return CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup(points, polygons) != 0;
}

template <typename PointRange, typename PolygonRange>
bool has_isolated_vertices(PointRange points, PolygonRange polygons)
{
  return CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(points, polygons) != 0;
}

template <typename PolygonMesh>
bool is_manifold(PolygonMesh& pmesh)
{
  return CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(pmesh, CGAL::parameters::dry_run(true)) == 0;
}

template <typename PolygonMesh>
bool is_self_intersecting(const PolygonMesh& pmesh) {
  return CGAL::Polygon_mesh_processing::does_self_intersect(pmesh);
}

template <typename PolygonMesh>
bool is_watertight(PolygonMesh& pmesh) {
  const bool manifold = is_manifold(pmesh);
  const bool closed = is_closed(pmesh);
  const bool self_intersecting = is_self_intersecting(pmesh);
  return manifold && closed && !self_intersecting;
}

template <typename PolygonMesh>
bool has_degenerate_faces(PolygonMesh& pmesh)
{
  std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor> dfs;
  CGAL::Polygon_mesh_processing::degenerate_faces(pmesh, std::inserter(dfs, dfs.begin()));
  return !dfs.empty();
}

template <typename PolygonMesh>
int connected_components(PolygonMesh& mesh)
{
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  auto fccmap = mesh.template add_property_map<face_descriptor, std::size_t>("f:CC").first;

  return CGAL::Polygon_mesh_processing::connected_components(mesh, fccmap);
}

template <typename PolygonMesh>
int euler_characteristic(const PolygonMesh& mesh)
{
  return mesh.number_of_vertices() - mesh.number_of_edges() + mesh.number_of_faces();
}

template <typename PolygonMesh>
int boundary_components(const PolygonMesh& mesh)
{
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  std::vector<halfedge_descriptor> border_cycles;

  CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));

  return border_cycles.size();
}

template <typename PolygonMesh>
int betti_0(PolygonMesh& mesh) {
    return connected_components(mesh);
}

template <typename PolygonMesh>
int betti_2(PolygonMesh& mesh) {
    return connected_components(mesh); // only for closed surfaces
}

template <typename PolygonMesh>
int betti_1(PolygonMesh& mesh) {
    return betti_0(mesh) + betti_2(mesh) - euler_characteristic(mesh);
}

// template <typename PolygonMesh>
// bool check_mesh_distance(const PolygonMesh& m0, const PolygonMesh& m1)
// {
//   auto dist = CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance<CGAL::Sequential_tag>(
//     m0, m1, CGAL::parameters::number_of_points_per_area_unit(4000));
//   std::cout << dist << std::endl;
//   return true;
// }

template<typename Domain>
void write_debug_grid(const Domain& domain, const std::string& filename) {
  using Point = typename Domain::Geom_traits::Point_3;
  using Mesh = CGAL::Surface_mesh<Point>;

  Mesh debug_grid;
  auto debug_grid_creator = [&](const typename Domain::cell_descriptor& c)
  {
    std::vector<typename Mesh::Vertex_index> cell_vertices;
    for (const auto& v : domain.cell_vertices(c)) {
      cell_vertices.push_back(debug_grid.add_vertex(domain.point(v)));
    }
    debug_grid.add_face(cell_vertices[6], cell_vertices[2], cell_vertices[0], cell_vertices[4]);
    debug_grid.add_face(cell_vertices[1], cell_vertices[3], cell_vertices[7], cell_vertices[5]);
    debug_grid.add_face(cell_vertices[0], cell_vertices[1], cell_vertices[5], cell_vertices[4]);
    debug_grid.add_face(cell_vertices[6], cell_vertices[7], cell_vertices[3], cell_vertices[2]);
    debug_grid.add_face(cell_vertices[2], cell_vertices[3], cell_vertices[1], cell_vertices[0]);
    debug_grid.add_face(cell_vertices[4], cell_vertices[5], cell_vertices[7], cell_vertices[6]);
  };
  domain.template for_each_cell<CGAL::Sequential_tag>(debug_grid_creator);
  CGAL::IO::write_OFF(filename, debug_grid);
}

#endif // CGAL_ISOSURFACING_TEST_UTIL_H
