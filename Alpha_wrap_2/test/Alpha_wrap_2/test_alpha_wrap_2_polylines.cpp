#define CGAL_AW2_TIMER
// #define CGAL_AW2_DEBUG_PP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Alpha_wrap_2/internal/validation.h>

#include <CGAL/bounding_box.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Random.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>

// #define CGAL_AW2_BREAK_ON_TEST_FAILURE

using namespace CGAL::Alpha_wraps_2::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = Kernel::FT;
using Point_2 = Kernel::Point_2;
using Vector_2 = Kernel::Vector_2;

using Points = std::vector<Point_2>;
using Polyline = std::vector<Point_2>;
using Polylines = std::vector<Polyline>;

using Multipolygon = CGAL::Multipolygon_with_holes_2<Kernel>;

CGAL::Bbox_2 compute_bbox(const Polylines& pls)
{
  CGAL::Bbox_2 bbox;
  for (const Polyline& pl : pls) {
    for (const Point_2& p : pl) {
      bbox += p.bbox();
    }
  }
  return bbox;
}

// Reads polylines from .wkt or .obj
bool read_polylines(const std::string& filename,
                    Polylines& polylines)
{
  std::ifstream in(filename);
  if(!in)
    return false;

  const std::string ext = CGAL::IO::internal::get_file_extension(filename);
  if(ext == "wkt")
  {
    // read_WKT() reads ALL multi-linestrings whereas read_multilinestring() reads only the first one
    Points pts;
    Multipolygon mp;
    return CGAL::IO::read_WKT(in, pts, polylines, mp);
  }
  else if (ext == "obj")
  {
    std::vector<Point_2> points;
    std::vector<std::vector<std::size_t> > id_polylines;
    std::vector<std::vector<std::size_t> > unused_id_polygons;
    bool success = CGAL::IO::internal::read_OBJ(in, points, id_polylines, unused_id_polygons);
    if(!success)
      return false;

    for(const std::vector<std::size_t>& id_pl : id_polylines) {
      polylines.emplace_back();
      for(const std::size_t pid : id_pl) {
        polylines.back().push_back(points[pid]);
      }
    }

    return !polylines.empty();
  }
  return false;
}

void read_translate_and_write(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Error: cannot open file '" << filename << "'" << std::endl;
    return;
  }

  Polylines polylines;
  bool res = read_polylines(filename, polylines);
  if(!res) {
    std::cerr << "Error: cannot read polylines from file '" << filename << "'" << std::endl;
    return;
  }

  std::cout << "Read " << polylines.size() << " polylines" << std::endl;

  std::string out_filename = filename;
  std::size_t dot_pos = out_filename.find_last_of('.');
  if (dot_pos != std::string::npos) {
    out_filename = out_filename.substr(0, dot_pos) + "_translated.wkt";
  } else {
    out_filename += "_translated.wkt";
  }

  // compute the bbox of all polylines
  CGAL::Bbox_2 bbox = compute_bbox(polylines);

  // Translate all points by 20% of the vertical size downwards
  // and by 50% of the horizontal size to the left
  const double translate_x = -0.5 * (bbox.xmax() - bbox.xmin());
  const double translate_y = -0.2 * (bbox.ymax() - bbox.ymin());
  const Vector_2 translation(translate_x, translate_y);

  Polylines translated_polylines = polylines;
  for(Polyline& pl : translated_polylines) {
    for(Point_2& p : pl) {
      p = p + translation;
    }
  }

  polylines.insert(polylines.end(), translated_polylines.begin(), translated_polylines.end());

  std::ofstream out(out_filename);
  if (!out) {
    std::cerr << "Error: cannot write to file '" << out_filename << "'" << std::endl;
    return;
  }

  out.precision(std::numeric_limits<double>::max_digits10);
  CGAL::IO::write_multi_linestring_WKT(out, polylines);

  std::cout << "Wrote " << polylines.size() << " translated polylines to '" << out_filename << "'" << std::endl;
}

// Wrap polylines, similar to alpha_wrap_point_set
bool alpha_wrap_polylines(Polylines& input_polylines,
                          const double alpha, const double offset)
{
  namespace AW2 = CGAL::Alpha_wraps_2;

  bool result = true;

  std::cout << "Input: " << input_polylines.size() << " polylines" << std::endl;
  std::cout << "Alpha: " << alpha << " Offset: " << offset << std::endl;

  const bool enforce_manifoldness = true;

  Multipolygon wrap;
  CGAL::alpha_wrap_2(input_polylines, alpha, offset, wrap,
                     CGAL::parameters::do_enforce_manifoldness(enforce_manifoldness));

  std::cout << "Result: " << wrap.polygons_with_holes().size() << " polygons" << std::endl;

  std::ofstream out_i("input.wkt");
  out_i.precision(std::numeric_limits<double>::max_digits10);
  CGAL::IO::write_multi_linestring_WKT(out_i, input_polylines);
  std::ofstream out_w("wrap.wkt");
  out_w.precision(std::numeric_limits<double>::max_digits10);
  CGAL::IO::write_multi_polygon_WKT(out_w, wrap);

  if(!AW2::internal::is_valid_wrap(wrap, enforce_manifoldness)) {
    std::cerr << "Error: invalid wrap" << std::endl;
    result = false;
  }

  if(!AW2::internal::is_outer_wrap_of_polylines(wrap, input_polylines)) {
    std::cerr << "Error: wrap does not contain the input" << std::endl;
    result = false;
  }

  if(!enforce_manifoldness) {
    if(!AW2::internal::has_bounded_edge_length(wrap, alpha)) {
      std::cerr << "Error: edge length check failure" << std::endl;
      result = false;
    }
  }

  return result;
}

// Compute alpha/offset and wrap polylines
bool alpha_wrap_polylines(const std::string& filename,
                          const double alpha_rel, const double offset_rel)
{
  Polylines input_polylines;
  bool res = read_polylines(filename, input_polylines);
  assert(res);
  assert(!input_polylines.empty());

  const CGAL::Bbox_2 bbox = compute_bbox(input_polylines);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()));

  const double alpha = diag_length / alpha_rel;
  const double offset = diag_length / offset_rel;

  return alpha_wrap_polylines(input_polylines, alpha, offset);
}

// Randomized test, similar to point version
bool alpha_wrap_polylines(const std::string& filename)
{
  bool result = true;

  CGAL::Random r;

  Polylines input_polylines;
  bool res = read_polylines(filename, input_polylines);
  assert(res);
  assert(!input_polylines.empty());

  const CGAL::Bbox_2 bbox = compute_bbox(input_polylines);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()));

  for(int i=0; i<3; ++i)
  {
    const double alpha_expo = r.get_double(3., 8);
    const double offset_expo = r.get_double(3., 8);
    const double alpha_rel = std::pow(2, alpha_expo);
    const double offset_rel = std::pow(2, offset_expo);
    const double alpha = diag_length / alpha_rel;
    const double offset = diag_length / offset_rel;

    std::cout << "===================================================" << std::endl;
    std::cout << filename << " " << alpha << " (rel " << alpha_rel << ")"
              << " " << offset << " (rel " << offset_rel << ")" << std::endl;
    std::cout << "Random seed = " << r.get_seed() << std::endl;

    if(!alpha_wrap_polylines(input_polylines, alpha, offset)) {
      result = false;
#ifdef CGAL_AW2_BREAK_ON_TEST_FAILURE
      break;
#endif
    }
  }

  return result;
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // read_translate_and_write(argv[1]);
  // return EXIT_SUCCESS;

  bool result = true;

  // For convenience to do manual testing
  if(argc > 1)
  {
    std::string arg1(argv[1]);
    if(std::filesystem::is_directory(arg1))
    {
      for(const auto& entry : std::filesystem::directory_iterator(arg1))
      {
        const std::string fname = entry.path().string();
        std::cout << "\n== " << fname << " ==" << std::endl;
        if(!alpha_wrap_polylines(fname)) {
          result = false;
#ifdef CGAL_AW2_BREAK_ON_TEST_FAILURE
          break;
#endif
        }
      }
    }
    else
    {
      const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
      const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

      result = alpha_wrap_polylines(arg1, relative_alpha, relative_offset);
    }
  }
  else
  {
    std::vector<std::string> default_inputs = { "data/polylines.wkt" };
    for(const std::string& filename : default_inputs) {
      if(!alpha_wrap_polylines(filename)) {
        result = false;
#ifdef CGAL_AW2_BREAK_ON_TEST_FAILURE
        break;
#endif
      }
    }
  }

  std::cout << "Done!" << std::endl;
  return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
