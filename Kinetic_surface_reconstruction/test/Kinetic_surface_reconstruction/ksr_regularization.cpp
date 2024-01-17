#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_surface_reconstruction_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/polygon_soup_io.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;
using Segment_3 = typename Kernel::Segment_3;

using Point_set = CGAL::Point_set_3<Point_3>;
using Point_map = typename Point_set::Point_map;
using Normal_map = typename Point_set::Vector_map;

using KSR = CGAL::Kinetic_surface_reconstruction_3<Kernel, Point_set, Point_map, Normal_map>;

int main(const int, const char**) {
  // Input.

  Point_set point_set;
  CGAL::IO::read_point_set("hilbert_cube.ply", point_set);

  auto with_reg = CGAL::parameters::maximum_distance(0.1)
    .maximum_angle(10)
    .k_neighbors(12)
    .minimum_region_size(10)
    .regularize_coplanarity(true)
    .regularize_parallelism(true)
    .maximum_offset(0.1)
    .angle_tolerance(10)
    .debug(true);

  auto without_reg = CGAL::parameters::maximum_distance(0.1)
    .maximum_angle(10)
    .k_neighbors(12)
    .minimum_region_size(10)
    .regularize_coplanarity(false)
    .regularize_parallelism(false);

  // Algorithm.
  KSR ksr(point_set);

  ksr.detect_planar_shapes(without_reg);

  std::size_t detected = ksr.detected_planar_shapes().size();

  ksr.detect_planar_shapes(with_reg);

  std::size_t regularized = ksr.detected_planar_shapes().size();

  std::cout << detected << " planar shapes regularized into " << regularized << std::endl;

  return regularized < detected ? EXIT_SUCCESS : EXIT_FAILURE;
}
