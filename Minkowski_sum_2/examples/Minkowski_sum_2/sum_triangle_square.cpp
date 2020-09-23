//! \file examples/Minkowski_sum_2/sum_triangle_square.cpp
// Computing the Minkowski sum of a triangle and a square.

#include <CGAL/minkowski_sum_2.h>

#include "bops_linear.h"

int main()
{
  // Construct the triangle.
  Polygon_2   P;
  P.push_back(Point_2(-1, -1));  P.push_back(Point_2(1, -1));
  P.push_back(Point_2(0, 1));
  std::cout << "P = " << P << std::endl;

  // Construct the square.
  Polygon_2   Q;
  Q.push_back(Point_2(3, -1));  Q.push_back(Point_2(5, -1));
  Q.push_back(Point_2(5, 1));   Q.push_back(Point_2(3,  1));
  std::cout << "Q = " << Q << std::endl;

  // Compute the Minkowski sum.
  Polygon_with_holes_2  sum = CGAL::minkowski_sum_2(P, Q);
  CGAL_assertion(sum.number_of_holes() == 0);
  std::cout << "P (+) Q = " << sum.outer_boundary() << std::endl;
  return 0;
}
