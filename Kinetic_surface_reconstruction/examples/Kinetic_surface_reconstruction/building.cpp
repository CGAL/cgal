#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_surface_reconstruction_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/PLY.h>
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
  std::string filename = CGAL::data_file_path("points_3/building.ply");

  Point_set point_set;
  std::ifstream input_file(filename, std::ios_base::binary);
  input_file >> point_set;
  input_file.close();

  auto param = CGAL::parameters::maximum_distance(0.1)
    .maximum_angle(10)
    .minimum_region_size(100)
    .reorient_bbox(true)
    .regularize_parallelism(true)
    .regularize_coplanarity(true)
    .angle_tolerance(5)
    .maximum_offset(0.02);

  // Algorithm.
  KSR ksr(point_set, param);

  ksr.detect_planar_shapes(param);

  ksr.initialize_partition(param);

  ksr.partition(2);

  std::vector<Point_3> vtx;
  std::vector<std::vector<std::size_t> > polylist;

  std::vector<FT> betas{0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99};

  for (FT b : betas) {
    vtx.clear();
    polylist.clear();

    ksr.reconstruct_with_ground(b, std::back_inserter(vtx), std::back_inserter(polylist));

    if (polylist.size() > 0)
      CGAL::IO::write_polygon_soup("polylist_" + std::to_string(b), vtx, polylist);
  }

  return EXIT_SUCCESS;
}
