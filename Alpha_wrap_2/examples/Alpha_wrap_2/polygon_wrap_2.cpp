#define CGAL_AW2_DEBUG_PP // @tmp (here and at other places)

#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

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

  const std::string filename = argc > 1 ? argv[1] : CGAL::data_file_path("wkt/issue.wkt");
  std::ifstream in(filename);
  Multipolygon mp_in;
  if(!in || !CGAL::IO::read_multi_polygon_WKT(in, mp_in))
  {
    std::cerr << "Can't read input file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << mp_in.polygons_with_holes().size() << " input polygons" << std::endl;

  // Optional:
  // use Polygon_repair to specify which strategy should be used to determine
  // what is inside and what is outside for invalid polygons.
  // We could also not repair, but then all edges (whether from the outer boundaries
  // or from the hole boundaries) are taken into account.
  auto rule = CGAL::Polygon_repair::Even_odd_rule();
  Multipolygon mp_repaired = CGAL::Polygon_repair::repair(mp_in, rule);

  std::cout << "post repair # = " << mp_repaired.polygons_with_holes().size() << std::endl;

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 10.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 300.;

  CGAL::Bbox_2 bbox = mp_repaired.bbox();
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()));
  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;
  std::cout << "absolute alpha = " << alpha << " absolute offset = " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  Multipolygon wrap;
  CGAL::alpha_wrap_2(mp_repaired, alpha, offset, wrap);

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
