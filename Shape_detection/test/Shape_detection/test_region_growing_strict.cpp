// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_segment_set.h>
#include <CGAL/Shape_detection/Region_growing/internal/free_functions.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_lines_points_with_normals() {

  using Point_2  = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;

  const std::vector< std::pair<Point_2, Vector_2> > points_with_normals = {
    std::make_pair(Point_2(0.1, 0.0), Vector_2(0.0, 1.0)),
    std::make_pair(Point_2(0.5, 0.0), Vector_2(0.0, 1.0)),
    std::make_pair(Point_2(0.9, 0.0), Vector_2(0.0, 1.0)),
    std::make_pair(Point_2(0.1, 2.0), Vector_2(0.0, 1.0)),
    std::make_pair(Point_2(0.5, 2.0), Vector_2(0.0, 1.0)),
    std::make_pair(Point_2(0.9, 2.0), Vector_2(0.0, 1.0)),
    std::make_pair(Point_2(0.1, 4.0), Vector_2(0.0, 1.0)),
    std::make_pair(Point_2(0.5, 4.0), Vector_2(0.0, 1.0)),
    std::make_pair(Point_2(0.9, 4.0), Vector_2(0.0, 1.0))
  };

  assert(points_with_normals.size() == 9);
  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_lines(
    points_with_normals, std::back_inserter(regions));
  assert(regions.size() == 3);
  assert(regions[0].size() == 3);
  assert(regions[1].size() == 3);
  assert(regions[2].size() == 3);
  return true;
}

template<class Kernel>
bool test_lines_polylines_2() {

  using Point_2 = typename Kernel::Point_2;
  const std::vector<Point_2> polyline_2 = {
    Point_2(0.10, 0.00), Point_2(0.50, 0.00), Point_2(0.90, 0.00),
    Point_2(0.13, 0.00), Point_2(0.17, 0.00), Point_2(0.21, 0.00),
    Point_2(0.21, 2.10), Point_2(0.21, 2.50), Point_2(0.21, 2.90),
    Point_2(0.21, 2.13), Point_2(0.21, 2.17), Point_2(0.21, 2.21)
  };

  assert(polyline_2.size() == 12);
  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_polylines(
    polyline_2, std::back_inserter(regions));
  assert(regions.size() == 2);
  assert(regions[0].size() == 6);
  assert(regions[1].size() == 6);
  return true;
}

template<class Kernel>
bool test_lines_polylines_3() {

  using Point_3 = typename Kernel::Point_3;
  const std::vector<Point_3> polyline_3 = {
    Point_3(0.10, 0.0, 1.0), Point_3(0.50, 0.0, 1.0), Point_3(0.90, 0.0, 1.0),
    Point_3(0.13, 0.0, 1.0), Point_3(0.17, 0.0, 1.0), Point_3(0.21, 0.0, 1.0)
  };

  assert(polyline_3.size() == 6);
  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_polylines(
    polyline_3, std::back_inserter(regions));
  assert(regions.size() == 1);
  assert(regions[0].size() == 6);
  return true;
}

template<class Kernel>
bool test_polylines_equal_points() {

  using Point_2 = typename Kernel::Point_2;
  const std::vector<Point_2> polyline_2 = {
    Point_2(0, 0), Point_2(1, 0), Point_2(2, 0), Point_2(3, 0),
    Point_2(7, 1), Point_2(8, 1), Point_2(9, 1), Point_2(10, 1),
    Point_2(14, 2), Point_2(15, 2), Point_2(16, 2), Point_2(17, 2),
    Point_2(19, 3), Point_2(20, 3), Point_2(21, 3), Point_2(22, 3)
  };

  assert(polyline_2.size() == 16);
  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_polylines(
    polyline_2, std::back_inserter(regions), CGAL::parameters::maximum_distance(0.01));
  assert(regions.size() == 4);

  assert(regions[0].size() == 4);
  assert(regions[1].size() == 4);
  assert(regions[2].size() == 4);
  assert(regions[3].size() == 4);
  return true;
}

template<class Kernel>
bool test_lines_segment_set_2() {

  using Point_2   = typename Kernel::Point_2;
  using Segment_2 = typename Kernel::Segment_2;

  using Segment_range = std::vector<Segment_2>;
  using Segment_map   = CGAL::Identity_property_map<Segment_2>;

  struct Neighbor_query {

    std::map<std::size_t, std::vector<std::size_t> > m_neighbors;
    Neighbor_query() {
      m_neighbors[0] = { 1 };
      m_neighbors[1] = {0, 2};
      m_neighbors[2] = {1, 3};
      m_neighbors[3] = { 2 };
    }

    void operator()(
      const std::size_t query_index, std::vector<std::size_t>& neighbors) const {
      neighbors.clear();
      const auto& data = m_neighbors.at(query_index);
      std::copy(data.begin(), data.end(), std::back_inserter(neighbors));
    }
  };

  using Region_type = CGAL::Shape_detection::
    Segment_set::Least_squares_line_fit_region<Kernel, Segment_range, Segment_map>;
  using Region_growing = CGAL::Shape_detection::
    Region_growing<Segment_range, Neighbor_query, Region_type>;

  const Segment_range segments = {
    Segment_2(Point_2(0.1, 0.0), Point_2(0.5, 0.0)),
    Segment_2(Point_2(0.5, 0.0), Point_2(0.9, 0.0)),
    Segment_2(Point_2(0.9, 0.0), Point_2(0.9, 0.5)),
    Segment_2(Point_2(0.9, 0.5), Point_2(0.9, 0.9))
  };
  assert(segments.size() == 4);

  Neighbor_query neighbor_query;
  Region_type region_type(segments);

  std::vector< std::vector<std::size_t> > regions;
  Region_growing region_growing(
    segments, neighbor_query, region_type);
  region_growing.detect(std::back_inserter(regions));
  assert(regions.size() == 2);
  assert(regions[0].size() == 2);
  assert(regions[1].size() == 2);
  return true;
}

template<class Kernel>
bool test_lines_segment_set_3() {

  using Point_3      = typename Kernel::Point_3;
  using Surface_mesh = CGAL::Surface_mesh<Point_3>;

  using Plane_region = CGAL::Shape_detection::
    Polygon_mesh::Least_squares_plane_fit_region<Kernel, Surface_mesh>;
  using Face_to_region_map = typename Plane_region::Face_to_region_map;

  using Polyline_graph = CGAL::Shape_detection::
    Polygon_mesh::Polyline_graph<Surface_mesh>;
  using Segment_range = typename Polyline_graph::Segment_range;
  using Segment_map = typename Polyline_graph::Segment_map;

  using Region_type = CGAL::Shape_detection::
    Segment_set::Least_squares_line_fit_region<Kernel, Segment_range, Segment_map>;
  using Sorting = CGAL::Shape_detection::
    Segment_set::Least_squares_line_fit_sorting<Kernel, Segment_range, Polyline_graph, Segment_map>;
  using Region_growing = CGAL::Shape_detection::
    Region_growing<Segment_range, Polyline_graph, Region_type, typename Sorting::Seed_map>;

  std::ifstream in(CGAL::data_file_path("meshes/am.off"));
  CGAL::set_ascii_mode(in);
  assert(in);

  Surface_mesh surface_mesh;
  in >> surface_mesh;
  in.close();

  assert(surface_mesh.number_of_faces() == 7320);
  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_planes(
    surface_mesh, std::back_inserter(regions));
  assert(regions.size() == 9);

  const auto face_range = faces(surface_mesh);
  assert(face_range.size() == 7320);

  const Face_to_region_map face_to_region_map(face_range, regions);
  Polyline_graph pgraph(surface_mesh, face_to_region_map);
  const auto& segment_range = pgraph.segment_range();

  Region_type region_type(
    segment_range, CGAL::parameters::segment_map(pgraph.segment_map()));
  Sorting sorting(
    segment_range, pgraph, CGAL::parameters::segment_map(pgraph.segment_map()));
  sorting.sort();

  regions.clear();
  Region_growing region_growing(
    segment_range, pgraph, region_type, sorting.seed_map());
  region_growing.detect(std::back_inserter(regions));
  assert(regions.size() == 21);
  return true;
}

template<class Kernel>
bool test_region_growing_lines() {

  assert(test_lines_points_with_normals<Kernel>());
  assert(test_lines_polylines_2<Kernel>());
  assert(test_lines_polylines_3<Kernel>());
  assert(test_polylines_equal_points<Kernel>());
  assert(test_lines_segment_set_2<Kernel>());
  assert(test_lines_segment_set_3<Kernel>());
  return true;
}

template<class Kernel>
bool test_planes_points_with_normals() {

  using Point_3  = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;

  const std::vector< std::pair<Point_3, Vector_3> > points_with_normals = {
    std::make_pair(Point_3(0.1, 0.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 0.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 0.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.1, 1.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 1.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 1.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.1, 2.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 2.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 2.0, 0.0), Vector_3(0.0, 0.0, 1.0))
  };

  assert(points_with_normals.size() == 9);
  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_planes(
    points_with_normals, std::back_inserter(regions));
  assert(regions.size() == 1);
  assert(regions[0].size() == 9);
  return true;
}

template<class Kernel>
bool test_planes_point_set() {

  using Point_3   = typename Kernel::Point_3;
  using Vector_3  = typename Kernel::Vector_3;
  using Point_set = CGAL::Point_set_3<Point_3, Vector_3>;

  std::vector< std::pair<Point_3, Vector_3> > points_with_normals = {
    std::make_pair(Point_3(0.1, 0.0, 1.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 0.0, 1.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 0.0, 1.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.1, 1.0, 1.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 1.0, 1.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 1.0, 1.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.1, 2.0, 1.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 2.0, 1.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 2.0, 1.0), Vector_3(0.0, 0.0, 1.0))
  };
  assert(points_with_normals.size() == 9);

  Point_set point_set(true);
  for (const auto& pwn : points_with_normals)
    point_set.insert(pwn.first, pwn.second);
  assert(point_set.size() == points_with_normals.size());

  points_with_normals.clear();
  assert(points_with_normals.size() == 0);

  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_planes(
    point_set, std::back_inserter(regions));
  assert(regions.size() == 1);
  assert(regions[0].size() == 9);
  return true;
}

template<class Kernel>
bool test_planes_polyhedron() {

  using Point_3    = typename Kernel::Point_3;
  using Polyhedron = CGAL::Polyhedron_3<
    Kernel, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>;

  Polyhedron polyhedron;
  const Point_3 p1(0, 0, 0);
  const Point_3 p2(1, 0, 0);
  const Point_3 p3(0, 1, 0);
  const Point_3 p4(0, 0, 1);
  const auto handle = polyhedron.make_tetrahedron(p1, p2, p3, p4);
  assert(polyhedron.is_tetrahedron(handle));

  assert(polyhedron.size_of_facets() == 4);
  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_planes(
    polyhedron, std::back_inserter(regions));
  assert(regions.size() == polyhedron.size_of_facets());
  return true;
}

template<class Kernel>
bool test_planes_surface_mesh() {

  using Point_3      = typename Kernel::Point_3;
  using Surface_mesh = CGAL::Surface_mesh<Point_3>;

  std::ifstream in(CGAL::data_file_path("meshes/am.off"));
  CGAL::set_ascii_mode(in);
  assert(in);

  Surface_mesh surface_mesh;
  in >> surface_mesh;
  in.close();

  assert(surface_mesh.number_of_faces() == 7320);
  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_planes(
    surface_mesh, std::back_inserter(regions));
  assert(regions.size() == 9);
  return true;
}

template<class Kernel>
bool test_region_growing_planes() {

  assert(test_planes_points_with_normals<Kernel>());
  assert(test_planes_point_set<Kernel>());
  assert(test_planes_polyhedron<Kernel>());
  assert(test_planes_surface_mesh<Kernel>());
  return true;
}

template<class Kernel>
bool test_region_growing_strict() {

  assert(test_region_growing_lines<Kernel>());
  assert(test_region_growing_planes<Kernel>());
  return true;
}

int main() {

  using SC = CGAL::Simple_cartesian<double>;
  using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
  using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

  // ------>

  bool sc_test_success = true;
  if (!test_region_growing_strict<SC>())
    sc_test_success = false;
  std::cout << "rg_strict, sc_test_success: " << sc_test_success << std::endl;
  assert(sc_test_success);

  // ------>

  bool epick_test_success = true;
  if (!test_region_growing_strict<EPICK>())
    epick_test_success = false;
  std::cout << "rg_strict, epick_test_success: " << epick_test_success << std::endl;
  assert(epick_test_success);

  // ------>

  bool epeck_test_success = true;
  if (!test_region_growing_strict<EPECK>())
    epeck_test_success = false;
  std::cout << "rg_strict, epeck_test_success: " << epeck_test_success << std::endl;
  assert(epeck_test_success);

  const bool success = sc_test_success && epick_test_success && epeck_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
