#ifndef CGAL_ISOSURFACING_TEST_UTIL_H
#define CGAL_ISOSURFACING_TEST_UTIL_H

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

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
bool has_degenerate_faces(PolygonMesh& pmesh)
{
  std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor> dfs;
  CGAL::Polygon_mesh_processing::degenerate_faces(pmesh, std::inserter(dfs, dfs.begin()));
  return !dfs.empty();
}

// template <typename PolygonMesh>
// bool check_mesh_distance(const PolygonMesh& m0, const PolygonMesh& m1)
// {
//   auto dist = CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance<CGAL::Sequential_tag>(
//     m0, m1, CGAL::parameters::number_of_points_per_area_unit(4000));
//   std::cout << dist << std::endl;
//   return true;
// }

#endif // CGAL_ISOSURFACING_TEST_UTIL_H
