#include <vector>
#include <cassert>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline.h>

// Typedefs.
using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Point_3  = typename Kernel::Point_3;

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/polyline_3.polylines.txt");
  CGAL::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Point_2> polyline_2;
  std::vector<Point_3> polyline_3;
  std::size_t n = std::size_t(-1);
  in >> n;
  polyline_2.reserve(n);
  polyline_3.reserve(n);
  while (n--) {
    FT a, b, c;
    in >> a >> b >> c;
    polyline_2.push_back(Point_2(a, b));
    polyline_3.push_back(Point_3(a, b, c));
    if (!in.good()) {
      // In case you have it, try to add en empty line at the end of the file.
      std::cout << "ERROR: cannot load the input polyline!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  in.close();

  std::cout << "* number of 2D input vertices: " << polyline_2.size() << std::endl;
  std::cout << "* number of 3D input vertices: " << polyline_3.size() << std::endl;
  assert(polyline_2.size() == 249);
  assert(polyline_3.size() == 249);

  // Run the algorithm both for 2D and 3D polyline.
  std::vector< std::vector<std::size_t> > segments_2;
  // CGAL::Shape_detection::region_growing_polylines(
  //   polyline_2, std::back_inserter(segments_2)
  //   CGAL::parameters::
  //   neighbor_radius(FT(5)).distance_threshold(FT(45) / FT(10)).
  //   angle_threshold(FT(45)).min_region_size(5));
  std::cout << "* number of found 2D segments: " << segments_2.size() << std::endl;
  assert(segments_2.size() == 65);

  std::vector< std::vector<std::size_t> > segments_3;
  // CGAL::Shape_detection::region_growing_polylines(
  //   polyline_3, std::back_inserter(segments_3)
  //   CGAL::parameters::
  //   neighbor_radius(FT(5)).distance_threshold(FT(45) / FT(10)).
  //   angle_threshold(FT(45)).min_region_size(5));
  std::cout << "* number of found 3D segments: " << segments_2.size() << std::endl;
  assert(segments_2.size() == 65);

  return EXIT_SUCCESS;
}
