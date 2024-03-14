#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_surface_reconstruction_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/polygon_soup_io.h>

using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT        = typename Kernel::FT;
using Point_3   = typename Kernel::Point_3;
using Vector_3  = typename Kernel::Vector_3;
using Segment_3 = typename Kernel::Segment_3;

using Point_set    = CGAL::Point_set_3<Point_3>;
using Point_map    = typename Point_set::Point_map;
using Normal_map   = typename Point_set::Vector_map;

using KSR = CGAL::Kinetic_surface_reconstruction_3<Kernel, Point_set, Point_map, Normal_map>;

int main(const int, const char**) {
  // Input.

  Point_set point_set;
  CGAL::IO::read_point_set("Cottage_cut.ply", point_set);

  auto param = CGAL::parameters::maximum_distance(0.33)
    .maximum_angle(32)
    .k_neighbors(12)
    .minimum_region_size(330)
    .reorient_bbox(false);

  // Algorithm.
  KSR ksr(point_set, param);

  ksr.detection_and_partition(1, param);

  std::vector<Point_3> vtx;
  std::vector<std::vector<std::size_t> > polylist;

  bool failed = false;

  ksr.reconstruct_with_ground(0.5, std::back_inserter(vtx), std::back_inserter(polylist));

  if (vtx.size() != 314 && polylist.size() != 91) {
    std::cerr << "reconstruction with reorientation and ground estimation set to outside provided wrong result: #vtx " << vtx.size() << " expected: 314, #polys: " << polylist.size() << " expected: 91" << std::endl;
    failed = true;
  }

  if (!failed)
    std::cout << "done!";

  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
