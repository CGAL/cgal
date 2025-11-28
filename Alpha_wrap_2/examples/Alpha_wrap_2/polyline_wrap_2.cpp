#define CGAL_AW2_DEBUG_PP // @tmp

#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>

#include <CGAL/Real_timer.h>
#include <CGAL/IO/WKT.h>

#include <iostream>
#include <string>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = K::Point_2;

using Points = std::vector<Point_2>;
using Polyline = std::vector<Point_2>;
using Polylines = std::vector<Polyline>;

using Multipolygon = CGAL::Multipolygon_with_holes_2<K>;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const std::string filename = argc > 1 ? argv[1] : CGAL::data_file_path("wkt/LetterAbis.wkt");

  // read_multi_linestring() expects an actual MULTILINESTRING entry whereas read_WKT() will read
  // all MULTILINESTRING and LINESTRING into a multi-linestring.
  Points pts_in;
  Polylines mls_in;
  Multipolygon mp_in;
  std::ifstream in(filename);
  if(!in || !CGAL::IO::read_WKT(in, pts_in, mls_in, mp_in))
  {
    std::cerr << "Can't read input file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << mls_in.size() << " input polylines" << std::endl;

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 300.;

  CGAL::Bbox_2 bbox;
  for(const Polyline& ls : mls_in) {
    for(const Point_2& pt : ls) {
      bbox += pt.bbox();
    }
  }

  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()));
  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;
  std::cout << "absolute alpha = " << alpha << " absolute offset = " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  Multipolygon wrap;
  CGAL::alpha_wrap_2(mls_in, alpha, offset, wrap);

  t.stop();
  std::cout << "Result: " << wrap.polygons_with_holes().size() << " polygons" << std::endl;
  std::cout << "Took " << t.time() << " s." << std::endl;

  // Save the result
  const std::string output_name = generate_output_name(filename, relative_alpha, relative_offset);
  std::cout << "Writing to " << output_name << std::endl;

  std::ofstream out(output_name);
  out.precision(std::numeric_limits<double>::max_digits10);
  CGAL::IO::write_multi_polygon_WKT(out, wrap);

  std::cout << "Done." << std::endl;

  return EXIT_SUCCESS;
}
