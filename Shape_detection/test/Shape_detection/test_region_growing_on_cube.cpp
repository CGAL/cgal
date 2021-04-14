// STL includes.
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/HalfedgeDS_vector.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>
#include <CGAL/Shape_detection/Region_growing/internal/free_functions.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_cube(int argc, char *argv[]) {

  using FT = typename Kernel::FT;

  using Polyhedron = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>;
  using Face_range = typename CGAL::Iterator_range<typename boost::graph_traits<Polyhedron>::face_iterator>;

  using Neighbor_query = SD::Polygon_mesh::One_ring_neighbor_query<Polyhedron, Face_range>;
  using Region_type    = SD::Polygon_mesh::Least_squares_plane_fit_region<Kernel, Polyhedron, Face_range>;
  using Sorting        = SD::Polygon_mesh::Least_squares_plane_fit_sorting<Kernel, Polyhedron, Neighbor_query, Face_range>;
  using Region_growing = SD::Region_growing<Face_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  using Vertex_to_point_map = typename Region_type::Vertex_to_point_map;

  // Default parameter values for the data file cube.off.
  const FT          distance_threshold = FT(1) / FT(10);
  const FT          angle_threshold    = FT(25);
  const std::size_t min_region_size    = 1;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/cube.off");
  CGAL::set_ascii_mode(in);
  assert(in);

  Polyhedron polyhedron;
  in >> polyhedron;
  in.close();
  const Face_range face_range = faces(polyhedron);
  assert(face_range.size() == 6);

  for (const auto& face : face_range) {
    const auto he = halfedge(face, polyhedron);
    const auto vertices = vertices_around_face(he, polyhedron);
    assert(vertices.size() == 4);
  }

  // Create parameter classes.
  Neighbor_query neighbor_query(polyhedron);
  const Vertex_to_point_map vertex_to_point_map(
    get(CGAL::vertex_point, polyhedron));

  Region_type region_type(
    polyhedron,
    CGAL::parameters::
    distance_threshold(distance_threshold).
    angle_threshold(angle_threshold).
    min_region_size(min_region_size).
    vertex_point_map(vertex_to_point_map));

  // Sort indices.
  Sorting sorting(
    polyhedron, neighbor_query,
    CGAL::parameters::vertex_point_map(vertex_to_point_map));
  sorting.sort();

  // Run region growing.
  Region_growing region_growing(
    face_range, neighbor_query, region_type, sorting.seed_map());

  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  assert(regions.size() == face_range.size());
  for (const auto& region : regions)
    assert(region_type.is_valid_region(region));

  std::vector<std::size_t> unassigned_faces;
  region_growing.unassigned_items(std::back_inserter(unassigned_faces));
  assert(unassigned_faces.size() == 0);

  // Test determenistic behavior and free functions.
  for (std::size_t k = 0; k < 3; ++k) {
    regions.clear();
    SD::internal::region_growing_planes(
      polyhedron, std::back_inserter(regions),
      CGAL::parameters::
      distance_threshold(distance_threshold).
      angle_threshold(angle_threshold).
      min_region_size(min_region_size).
      vertex_point_map(vertex_to_point_map));
    assert(regions.size() == face_range.size());
  }
  return true;
}

int main(int argc, char *argv[]) {

  using SC = CGAL::Simple_cartesian<double>;
  using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
  using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

  // ------>

  bool sc_test_success = true;
  if (!test_region_growing_on_cube<SC>(argc, argv))
    sc_test_success = false;
  std::cout << "rg_cube, sc_test_success: " << sc_test_success << std::endl;
  assert(sc_test_success);

  // ------>

  bool epick_test_success = true;
  if (!test_region_growing_on_cube<EPICK>(argc, argv))
    epick_test_success = false;
  std::cout << "rg_cube, epick_test_success: " << epick_test_success << std::endl;
  assert(epick_test_success);

  // ------>

  bool epeck_test_success = true;
  if (!test_region_growing_on_cube<EPECK>(argc, argv))
    epeck_test_success = false;
  std::cout << "rg_cube, epeck_test_success: " << epeck_test_success << std::endl;
  assert(epeck_test_success);

  const bool success = sc_test_success && epick_test_success && epeck_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
