// STL includes.
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Polygon_mesh.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_polygon_mesh(int argc, char *argv[]) {

  using FT      = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
  using Face_range   = typename Polygon_mesh::Face_range;

  using Neighbor_query = SD::Polygon_mesh::One_ring_neighbor_query<Polygon_mesh>;
  using Region_type    = SD::Polygon_mesh::Least_squares_plane_fit_region<Kernel, Polygon_mesh>;
  using Region_growing = SD::Region_growing<Neighbor_query, Region_type>;

  // Default parameter values.
  const FT          distance_threshold = FT(1);
  const FT          angle_threshold    = FT(45);
  const std::size_t min_region_size    = 1;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : CGAL::data_file_path("meshes/building.off"));
  CGAL::IO::set_ascii_mode(in);
  assert(in);

  Polygon_mesh mesh;
  in >> mesh;
  in.close();
  const Face_range face_range = faces(mesh);
  assert(face_range.size() == 32245);

  using Vertex_to_point_map = typename Region_type::Vertex_to_point_map;
  const Vertex_to_point_map vertex_to_point_map(
    get(CGAL::vertex_point, mesh));

  // Create parameter classes.
  Neighbor_query neighbor_query(mesh);
  Region_type region_type(
    mesh,
    CGAL::parameters::
    maximum_distance(distance_threshold).
    maximum_angle(angle_threshold).
    minimum_region_size(min_region_size).
    vertex_point_map(vertex_to_point_map));

  // Run region growing.
  Region_growing region_growing(
    face_range, neighbor_query, region_type);

  std::vector<typename Region_growing::Primitive_and_region> regions;
  region_growing.detect(std::back_inserter(regions));
  assert(regions.size() == 1077);
  for (const auto& region : regions)
    assert(region_type.is_valid_region(region.second));

  auto map = region_growing.region_map();
  for (auto fit : face_range) {
    std::size_t id = get(region_growing.region_map(), fit);
    assert(id != std::size_t(-1));
  }

  std::vector<typename Region_growing::Item> unassigned_faces;
  region_growing.unassigned_items(face_range, std::back_inserter(unassigned_faces));
  assert(unassigned_faces.size() == 0);
  return true;
}

int main(int argc, char *argv[]) {

  using SC = CGAL::Simple_cartesian<double>;
  using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
  using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

  // ------>

  bool sc_test_success = true;
  if (!test_region_growing_on_polygon_mesh<SC>(argc, argv))
    sc_test_success = false;
  std::cout << "rg_pmesh, sc_test_success: " << sc_test_success << std::endl;
  assert(sc_test_success);

  // ------>

  bool epick_test_success = true;
  if (!test_region_growing_on_polygon_mesh<EPICK>(argc, argv))
    epick_test_success = false;
  std::cout << "rg_pmesh, epick_test_success: " << epick_test_success << std::endl;
  assert(epick_test_success);

  // ------>

  bool epeck_test_success = true;
  if (!test_region_growing_on_polygon_mesh<EPECK>(argc, argv))
    epeck_test_success = false;
  std::cout << "rg_pmesh, epeck_test_success: " << epeck_test_success << std::endl;
  assert(epeck_test_success);

  const bool success = sc_test_success && epick_test_success && epeck_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
