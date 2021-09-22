// Turn off warnings.
#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

// STL includes.
#include <string>
#include <vector>
#include <iostream>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>

namespace SD = CGAL::Shape_detection;
using Kernel = CGAL::Simple_cartesian<double>;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Point_3  = typename Kernel::Point_3;
using Vector_2 = typename Kernel::Vector_2;
using Vector_3 = typename Kernel::Vector_3;

using Point_with_normal_2 = std::pair<Point_2, Vector_2>;
using Point_with_normal_3 = std::pair<Point_3, Vector_3>;

using Input_range_2 = std::vector<Point_with_normal_2>;
using Input_range_3 = std::vector<Point_with_normal_3>;

using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
using Face_range   = typename Polygon_mesh::Face_range;

using Point_map_2 = CGAL::First_of_pair_property_map<Point_with_normal_2>;
using Point_map_3 = CGAL::First_of_pair_property_map<Point_with_normal_3>;

using Normal_map_2 = CGAL::Second_of_pair_property_map<Point_with_normal_2>;
using Normal_map_3 = CGAL::Second_of_pair_property_map<Point_with_normal_3>;

using NQ1 = SD::Point_set::K_neighbor_query<Kernel, Input_range_2, Point_map_2>;
using NQ2 = SD::Point_set::Sphere_neighbor_query<Kernel, Input_range_3, Point_map_3>;
using NQ3 = SD::Polygon_mesh::One_ring_neighbor_query<Polygon_mesh>;

using RT1 = SD::Point_set::Least_squares_line_fit_region<Kernel, Input_range_2, Point_map_2, Normal_map_2>;
using RT2 = SD::Point_set::Least_squares_plane_fit_region<Kernel, Input_range_3, Point_map_3, Normal_map_3>;
using RT3 = SD::Polygon_mesh::Least_squares_plane_fit_region<Kernel, Polygon_mesh>;

using S1 = SD::Point_set::Least_squares_line_fit_sorting<Kernel, Input_range_2, NQ1, Point_map_2>;
using S2 = SD::Point_set::Least_squares_plane_fit_sorting<Kernel, Input_range_3, NQ2, Point_map_3>;
using S3 = SD::Polygon_mesh::Least_squares_plane_fit_sorting<Kernel, Polygon_mesh, NQ3>;

using RG1 = SD::Region_growing<Input_range_2, NQ1, RT1, typename S1::Seed_map>;
using RG2 = SD::Region_growing<Input_range_3, NQ2, RT2, typename S2::Seed_map>;
using RG3 = SD::Region_growing<Face_range   , NQ3, RT3, typename S3::Seed_map>;

int main() {

  // Default parameter values.
  const std::size_t k                  = 6;
  const FT          distance_threshold = FT(2);
  const FT          angle_threshold    = FT(15);
  const std::size_t min_region_size    = 1;
  const FT          sphere_radius      = FT(2);

  const Kernel traits;
  const Input_range_2 input_range_2 = {
    std::make_pair(Point_2(0, 0), Vector_2(0, 1)),
    std::make_pair(Point_2(1, 0), Vector_2(0, 1))
  };
  const Input_range_3 input_range_3 = {
    std::make_pair(Point_3(0, 0, 0), Vector_3(0, 0, 1)),
    std::make_pair(Point_3(1, 0, 0), Vector_3(0, 0, 1)),
    std::make_pair(Point_3(2, 0, 0), Vector_3(0, 0, 1))
  };
  Polygon_mesh polygon_mesh;
  CGAL::make_tetrahedron(
    Point_3(0, 0, 0), Point_3(2, 0, 0),
    Point_3(1, 1, 1), Point_3(1, 0, 2), polygon_mesh);
  const auto vertex_to_point_map =
    get_const_property_map(CGAL::vertex_point, polygon_mesh);
  const auto face_range = faces(polygon_mesh);

  Point_map_2 point_map_2;
  Point_map_3 point_map_3;
  Normal_map_2 normal_map_2;
  Normal_map_3 normal_map_3;

  assert(input_range_2.size() == 2);
  assert(input_range_3.size() == 3);
  assert(polygon_mesh.number_of_faces() == 4);

  NQ1 nq1(input_range_2, k, point_map_2);
  NQ2 nq2(input_range_3, sphere_radius, point_map_3);
  NQ3 nq3(polygon_mesh);

  RT1 rt1(input_range_2,
    distance_threshold, angle_threshold, min_region_size,
    point_map_2, normal_map_2, traits);
  RT2 rt2(input_range_3,
    distance_threshold, angle_threshold, min_region_size,
    point_map_3, normal_map_3, traits);
  RT3 rt3(polygon_mesh,
    distance_threshold, angle_threshold, min_region_size,
    vertex_to_point_map, traits);

  S1 s1(input_range_2, nq1, point_map_2);
  S2 s2(input_range_3, nq2, point_map_3);
  S3 s3(polygon_mesh , nq3, vertex_to_point_map);

  s1.sort();
  s2.sort();
  s3.sort();

  RG1 rg1(input_range_2, nq1, rt1, s1.seed_map());
  RG2 rg2(input_range_3, nq2, rt2, s2.seed_map());
  RG3 rg3(face_range   , nq3, rt3, s3.seed_map());

  std::vector< std::vector<std::size_t> > regions1;
  std::vector< std::vector<std::size_t> > regions2;
  std::vector< std::vector<std::size_t> > regions3;
  rg1.detect(std::back_inserter(regions1));
  rg2.detect(std::back_inserter(regions2));
  rg3.detect(std::back_inserter(regions3));

  assert(regions1.size() == 1);
  assert(regions2.size() == 1);
  assert(regions3.size() == 4);

  std::cout << "rg_deprecated, sc_test_success: " << true << std::endl;
  return EXIT_SUCCESS;
}
