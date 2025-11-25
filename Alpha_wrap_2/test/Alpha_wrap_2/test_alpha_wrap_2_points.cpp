#define CGAL_AW3_TIMER

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Alpha_wrap_2/internal/validation.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/bounding_box.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Random.h>

#include <fstream>
#include <iostream>
#include <vector>

using namespace CGAL::Alpha_wraps_2::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = Kernel::FT;
using Point_2 = Kernel::Point_2;
using Vector_2 = Kernel::Vector_2;
using Iso_rectangle_2 = Kernel::Iso_rectangle_2;

using Points = std::vector<Point_2>;

using Multipolygon = CGAL::Multipolygon_with_holes_2<Kernel>;

// small utility to construct a multipoint .wkt from a .xy or .xyz
void xyz_to_multipoint(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Error: cannot open file '" << filename << "'" << std::endl;
    return;
  }

  Points points;
  std::string line;
  int lineNumber = 0;
  while (std::getline(in, line)) {
    ++lineNumber;
    // Trim whitespace
    line.erase(line.find_last_not_of(" \t\r\n") + 1);
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    // Skip empty or comment lines
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream iss(line);
    std::vector<double> values;
    double x, y;
    // ignore everything except the (X,Y) coordinates
    if (iss >> x >> y) {
      points.emplace_back(x, y);
    }
  }

  std::cout << "Read " << points.size() << " points" << std::endl;

  std::string out_filename = filename;
  std::size_t dot_pos = out_filename.find_last_of('.');
  if (dot_pos != std::string::npos) {
    out_filename = out_filename.substr(0, dot_pos) + ".wkt";
  } else {
    out_filename += ".wkt";
  }

  std::ofstream out(out_filename);
  if (!out) {
    std::cerr << "Error: cannot write to file '" << out_filename << "'" << std::endl;
    return;
  }
  out.precision(std::numeric_limits<double>::max_digits10);
  CGAL::IO::write_multi_point_WKT(out, points);

  std::cout << "Wrote " << points.size() << " points to '" << out_filename << "'" << std::endl;
}

void alpha_wrap_point_set(Points& input_points,
                          const double alpha,
                          const double offset)
{
  namespace AW2 = CGAL::Alpha_wraps_2;
  namespace PMP = CGAL::Polygon_mesh_processing;

  std::cout << "Input: " << input_points.size() << " points" << std::endl;
  std::cout << "Alpha: " << alpha << " Offset: " << offset << std::endl;

  const bool enforce_manifoldness = true;

  Multipolygon wrap;
  CGAL::alpha_wrap_2(input_points, alpha, offset, wrap,
                     CGAL::parameters::do_enforce_manifoldness(enforce_manifoldness));

  std::cout << "Result: " << wrap.polygons_with_holes().size() << " polygons" << std::endl;

  std::ofstream out("last.off");
  out.precision(17);
  CGAL::IO::write_multi_polygon_WKT(out, wrap);

  assert(AW2::internal::is_valid_wrap(wrap, enforce_manifoldness));
  assert(AW2::internal::is_outer_wrap_of_point_set(wrap, input_points));

  // if(!enforce_manifoldness)
  //   assert(AW2::internal::has_expected_Hausdorff_distance(wrap, input_points, alpha, offset));

  if(!enforce_manifoldness)
    assert(AW2::internal::has_bounded_edge_length(wrap, alpha));
}

void alpha_wrap_point_set(const std::string& filename,
                          const double alpha_rel,
                          const double offset_rel)
{
  std::ifstream in(filename);
  Points input_points;
  bool res = CGAL::IO::read_multi_point_WKT(in, input_points);
  assert(res);
  assert(!input_points.empty());

  const Iso_rectangle_2 bbox = CGAL::bounding_box(input_points.begin(), input_points.end());
  const Vector_2 longest_diag = Point_2(bbox.xmax(), bbox.ymax()) -
                                Point_2(bbox.xmin(), bbox.ymin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  const double alpha = longest_diag_length / alpha_rel;
  const double offset = longest_diag_length / offset_rel;

  alpha_wrap_point_set(input_points, alpha, offset);
}

void alpha_wrap_point_set(const std::string& filename)
{
  std::cout << "\n== " << filename << "==" << std::endl;

  CGAL::Random r;

  std::ifstream in(filename);
  assert(in);
  Points input_points;
  bool res = CGAL::IO::read_multi_point_WKT(in, input_points);
  assert(res);
  assert(!input_points.empty());

  const Iso_rectangle_2 bbox = CGAL::bounding_box(input_points.begin(), input_points.end());
  const Vector_2 longest_diag = Point_2(bbox.xmax(), bbox.ymax()) -
                                Point_2(bbox.xmin(), bbox.ymin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  for(int i=0; i<2; ++i)
  {
    const double alpha_expo = r.get_double(0., 6); // to have alpha_rel between 1 and 64
    const double offset_expo = r.get_double(0., 6);
    const double alpha_rel = std::pow(2, alpha_expo);
    const double offset_rel = std::pow(2, offset_expo);
    const double alpha = longest_diag_length / alpha_rel;
    const double offset = longest_diag_length / offset_rel;

    std::cout << "===================================================" << std::endl;
    std::cout << filename << " " << alpha << " (rel " << alpha_rel << ")"
                          << " " << offset << " (rel " << offset_rel << ")" << std::endl;
    std::cout << "Random seed = " << r.get_seed() << std::endl;

    alpha_wrap_point_set(input_points, alpha, offset);
  }
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // xyz_to_multipoint(argv[1]);
  // return EXIT_SUCCESS;

  // For convenience to do manual testing
  if(argc > 1)
  {
    const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
    const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

    alpha_wrap_point_set(argv[1], relative_alpha, relative_offset);

    std::cout << "Done!" << std::endl;
    return EXIT_SUCCESS;
  }

  alpha_wrap_point_set("data/buildings_outline.wkt");
  alpha_wrap_point_set("data/circles.wkt");
  alpha_wrap_point_set("data/duplicate_points.wkt");
  alpha_wrap_point_set("data/four_points.wkt");
  alpha_wrap_point_set("data/hippo2.wkt");
  alpha_wrap_point_set("data/kitten.wkt");
  alpha_wrap_point_set("data/spheres.wkt");

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
