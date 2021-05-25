#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

using Kernel  = CGAL::Simple_cartesian<double>;
using Point_2 = typename Kernel::Point_2;

int main() {

  // Create input contour.
  const std::vector<Point_2> contour = {
    Point_2(0.00,  0.00),
    Point_2(0.50, -0.05),
    Point_2(1.00,  0.00),
    Point_2(1.05,  0.50),
    Point_2(1.00,  1.00),
    Point_2(0.00,  1.00)
  };

  // Regularize this contour.
  std::vector<Point_2> regularized;
  CGAL::Shape_regularization::Contours::
    regularize_closed_contour(contour, std::back_inserter(regularized));
}
