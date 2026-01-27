#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_repair/repair.h>

#include <CGAL/Real_timer.h>
#include <CGAL/IO/WKT.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = K::Point_2;

using Multipolygon = CGAL::Multipolygon_with_holes_2<K>;
using Polygon_with_holes = Multipolygon::Polygon_with_holes_2;
using Polygon = Polygon_with_holes::Polygon_2;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // Read the input multipolygon
  const std::string filename = argc > 1 ? argv[1] : CGAL::data_file_path("wkt/issue.wkt");
  std::ifstream in(filename);
  Multipolygon mp_in;
  if(!in || !CGAL::IO::read_multi_polygon_WKT(in, mp_in))
  {
    std::cerr << "Can't read input file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << mp_in.polygons_with_holes().size() << " input polygons" << std::endl;

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 10.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 300.;

  CGAL::Bbox_2 bbox = mp_in.bbox();
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()));
  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;
  std::cout << "absolute alpha = " << alpha << " absolute offset = " << offset << std::endl;

  // Set up seed points
  //
  // There is no limit on how many seeds can be used.
  // However, the algorithm automatically determines whether a seed can be used
  // to initialize the refinement based on a few conditions (distance to the offset, value of alpha, etc.)
  // See internal function Alpha_wrap_3::initialize_from_cavities() for more information
  std::vector<Point_2> seeds = {
    // the center of the bbox
    Point_2((bbox.xmin() + bbox.xmax()) / 2.0, (bbox.ymin() + bbox.ymax()) / 2.0),

    // a point outside the bounding box
    Point_2(bbox.xmax() + 10 * (alpha + offset),
            bbox.ymax() + 10 * (alpha + offset))
  };

  // Construct the wrap using the seed point
  CGAL::Real_timer t;
  t.start();

  Multipolygon wrap;
  CGAL::alpha_wrap_2(mp_in, alpha, offset, wrap, CGAL::parameters::seed_points(std::ref(seeds)));

  t.stop();
  std::cout << "Result: " << wrap.polygons_with_holes().size() << " polygons" << std::endl;
  std::cout << "Took " << t.time() << " s." << std::endl;

  // Save the result
  std::ofstream out("cavity_wrap.wkt");
  out.precision(std::numeric_limits<double>::max_digits10);
  CGAL::IO::write_multi_polygon_WKT(out, wrap);

  std::cout << "Done." << std::endl;

  return EXIT_SUCCESS;
}
