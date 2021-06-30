#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

// Typedefs.
using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Contour = std::vector<Point_2>;

using Contour_directions =
  CGAL::Shape_regularization::Contours::Multiple_directions_2<Kernel, Contour>;

int main(int argc, char *argv[]) {

  // If we want to load a different file, we load it from a path.
  std::string path = "data/contour.polylines";
  if (argc > 1) path = argv[1];
  Saver<Kernel> saver;

  // Set parameters.
  const FT min_length_2 = FT(2);
  const FT  max_angle_2 = FT(20);
  const FT max_offset_2 = FT(2);

  // Initialize contour.
  std::vector<Point_2> contour;
  initialize_contour(path, contour);

  // Save input contour.
  if (argc > 2) {
    const std::string full_path = std::string(argv[2]) + "regularize_closed_contour_before";
    saver.export_eps_closed_contour(contour, full_path, FT(8));
  }

  // Regularize.
  const bool is_closed = true;
  Contour_directions directions(
    contour, is_closed, CGAL::parameters::
    minimum_length(min_length_2).maximum_angle(max_angle_2));

  std::vector<Point_2> regularized;
  CGAL::Shape_regularization::Contours::regularize_closed_contour(
    contour, directions, std::back_inserter(regularized),
    CGAL::parameters::maximum_offset(max_offset_2));

  std::cout << "* number of directions = " <<
    directions.number_of_directions() << std::endl;

  // Save regularized contour.
  if (argc > 2) {
    const std::string full_path = std::string(argv[2]) + "regularize_closed_contour_after";
    saver.export_eps_closed_contour(regularized, full_path, FT(8));
  }
}
