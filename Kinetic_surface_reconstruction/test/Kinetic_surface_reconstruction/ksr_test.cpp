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

int main(const int argc, const char** argv) {
  // Input.

  Point_set point_set;
  CGAL::IO::read_point_set(CGAL::data_file_path("points_3/building.xyz"), point_set);

  auto param = CGAL::parameters::maximum_distance(0.05)
    .maximum_angle(10)
    .k_neighbors(12)
    .minimum_region_size(50);

  // Algorithm.
  KSR ksr(point_set, param);

  ksr.detection_and_partition(1, param);

  std::vector<Point_3> vtx;
  std::vector<std::vector<std::size_t> > polylist;

  std::map<typename KSR::KSP::Face_support, bool> external_nodes;

  bool failed = false;

  ksr.reconstruct(0.5, external_nodes, std::back_inserter(vtx), std::back_inserter(polylist));

  if (vtx.size() != 10 && polylist.size() != 7) {
    std::cerr << "reconstruction with external nodes set to outside provided wrong result: #vtx " << vtx.size() << " expected: 10, #polys: " << polylist.size() << " expected: 7" << std::endl;
    failed = true;
  }

  vtx.clear();
  polylist.clear();

  ksr.reconstruct_with_ground(0.5, std::back_inserter(vtx), std::back_inserter(polylist));

  if (vtx.size() != 22 && polylist.size() != 7) {
    std::cerr << "reconstruction with ground set to outside provided wrong result: #vtx " << vtx.size() << " expected: 22, #polys: " << polylist.size() << " expected: 7" << std::endl;
    failed = true;
  }

  vtx.clear();
  polylist.clear();

  ksr.reconstruct_with_ground(0.999, std::back_inserter(vtx), std::back_inserter(polylist));

  if (vtx.size() != 0 && polylist.size() != 0) {
    std::cerr << "reconstruction with high lambda provided wrong result: #vtx " << vtx.size() << " expected: 0, #polys: " << polylist.size() << " expected: 0" << std::endl;
    failed = true;
  }

  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
