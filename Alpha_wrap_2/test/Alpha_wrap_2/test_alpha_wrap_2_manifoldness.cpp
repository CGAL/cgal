#define CGAL_AW2_TIMER
#define CGAL_AW2_DEBUG
#define CGAL_AW2_DEBUG_MANIFOLDNESS
//#define CGAL_AW2_DEBUG_STEINER_COMPUTATION
//#define CGAL_AW2_DEBUG_INITIALIZATION
//#define CGAL_AW2_DEBUG_QUEUE

#include "test_utilities.h"

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Alpha_wrap_2/internal/validation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/draw_multipolygon_with_holes_2.h>
#include <CGAL/Random.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>

using namespace CGAL::Alpha_wraps_2::internal;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = K::FT;
using Point_2 = K::Point_2;
using Vector_2 = K::Vector_2;

using Points = std::vector<Point_2>;
using Polyline = std::vector<Point_2>;
using Polylines = std::vector<Polyline>;
using Multipolygon = CGAL::Multipolygon_with_holes_2<K>;

bool alpha_wrap_triangle_manifoldness(Polylines& input_polylines,
                                      const double alpha,
                                      const double offset)
{
  namespace AW2 = CGAL::Alpha_wraps_2;

  std::cout << "Input: " << input_polylines.size() << " polylines" << std::endl;
  std::cout << "Alpha: " << alpha << " Offset: " << offset << std::endl;

  std::cout << "  --- WITHOUT Manifoldness Enforcement" << std::endl;

  Multipolygon nm_wrap;
  CGAL::alpha_wrap_2(input_polylines, alpha, offset, nm_wrap,
                     CGAL::parameters::do_enforce_manifoldness(false));
  CGAL::draw(nm_wrap);

  assert(AW2::internal::is_valid_wrap(nm_wrap, false /*manifoldness*/));
  assert(AW2::internal::is_outer_wrap_of_polylines(nm_wrap, input_polylines));
  assert(AW2::internal::has_bounded_edge_length(nm_wrap, alpha));
  // assert(AW2::internal::has_expected_Hausdorff_distance(nm_wrap, input_polylines, alpha, offset));

  FT base_area = area(nm_wrap);

  std::cout << "  --- WITH Manifoldness Enforcement" << std::endl;

  Multipolygon m_wrap;
  CGAL::alpha_wrap_2(input_polylines, alpha, offset, m_wrap,
                     CGAL::parameters::do_enforce_manifoldness(true));
  CGAL::draw(m_wrap);

  assert(AW2::internal::is_valid_wrap(m_wrap, true /*manifoldness*/));
  assert(AW2::internal::is_outer_wrap_of_polylines(m_wrap, input_polylines));
  // enforcing manifoldness may slightly change the Hausdorff distance,
  // so this can't be enforced

  const FT final_area = area(m_wrap);

  if(base_area != 0)
  {
    const FT ratio = final_area / base_area;

    std::cout << "Areas post-manifoldness fix:\n"
              << "before: " << base_area << "\n"
              << "after:  " << final_area << "\n"
              << "ratio:  " << ratio << std::endl;
    if(ratio > 1.1) { // more than 10% extra area
      std::cerr << "W: large increase of area after manifoldness resolution" << std::endl;
      return false;
    }
  }

  return true;
}

bool alpha_wrap_triangle_manifoldness(const std::string& filename,
                                      const double alpha_rel,
                                      const double offset_rel)
{
  std::cout << "\n===================================================" << std::endl;
  std::cout << "FILE: " << filename << std::endl;

  Polylines input_polylines;
  bool res = read_polylines(filename, input_polylines);
  assert(res);
  assert(!input_polylines.empty());

  CGAL::Bbox_2 bbox = compute_bbox(input_polylines);
  const Vector_2 longest_diag = Point_2(bbox.xmax(), bbox.ymax()) -
                                Point_2(bbox.xmin(), bbox.ymin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  const double alpha = longest_diag_length / alpha_rel;
  const double offset = longest_diag_length / offset_rel;

  return alpha_wrap_triangle_manifoldness(input_polylines, alpha, offset);
}

bool alpha_wrap_triangle_manifoldness(const std::string& filename)
{
  std::cout << "\n===================================================" << std::endl;
  std::cout << "FILE: " << filename << std::endl;

  CGAL::Random r;
  std::cout << "Random seed = " << r.get_seed() << std::endl;

  Polylines input_polylines;
  bool res = read_polylines(filename, input_polylines);
  assert(res);
  assert(!input_polylines.empty());

  CGAL::Bbox_2 bbox = compute_bbox(input_polylines);
  const Vector_2 longest_diag = Point_2(bbox.xmax(), bbox.ymax()) -
                                Point_2(bbox.xmin(), bbox.ymin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  const double alpha_expo = r.get_double(0., 6); // to have alpha_rel between 1 and 64
  const double offset_expo = r.get_double(0., 6);
  const double alpha_rel = std::pow(2, alpha_expo);
  const double offset_rel = std::pow(2, offset_expo);
  const double alpha = longest_diag_length / alpha_rel;
  const double offset = longest_diag_length / offset_rel;

  return alpha_wrap_triangle_manifoldness(input_polylines, alpha, offset);
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  bool result = true;

  if(argc > 1)
  {
    std::string arg1(argv[1]);
    if(std::filesystem::is_directory(arg1))
    {
      for(const auto& entry : std::filesystem::directory_iterator(arg1))
      {
        const std::string fname = entry.path().string();
        std::cout << "\n== " << fname << " ==" << std::endl;

        bool local_result;
        if (argc > 2) {
          const double relative_alpha = std::stod(argv[2]);
          const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;
          local_result = alpha_wrap_triangle_manifoldness(fname, relative_alpha, relative_offset);
        } else {
          local_result = alpha_wrap_triangle_manifoldness(fname);
        }

        if(!local_result) {
          result = false;
          std::cerr << "Error!" << std::endl;
#ifdef CGAL_AW2_BREAK_ON_TEST_FAILURE
          break;
#endif
        }
        std::cout << "OK" << std::endl;
      }
    }
    else
    {
      const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
      const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

      result = alpha_wrap_triangle_manifoldness(arg1, relative_alpha, relative_offset);
    }
  }
  else
  {
    std::vector<std::tuple<std::string, double, double> > default_inputs =
      {
        {"data/stories.obj", 5, 600},
        {"data/heart.obj", 7, 600},
        {"data/halloween_006.obj", 10, 600},
        {"data/duel.obj", 12, 600},
        {"data/Letra.obj", 15, 600},
        {"data/butterflies.obj", 18, 600},
        {"data/Easter.obj", 20, 600},
      };

    for(const auto& entry : default_inputs) {
      if(!alpha_wrap_triangle_manifoldness(std::get<0>(entry),
                                           std::get<1>(entry),
                                           std::get<2>(entry))) {
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
