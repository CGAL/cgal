//! \file examples/Minkowski_sum_2/sum_triangle_square.cpp
// Computing the Minkowski sum of a triangle and a square.

#include <CGAL/basic.h>
#include <CGAL/minkowski_sum_2.h>

#include "bops_linear.h"

int main()
{
  // Construct the triangle.
  Polygon   P;
  P.push_back(Point(-1, -1));  P.push_back(Point(1, -1));
  P.push_back(Point(0, 1));
  std::cout << "P = " << P << std::endl;

  // Construct the square.
  Polygon   Q;
  Q.push_back(Point(3, -1));  Q.push_back(Point(5, -1));
  Q.push_back(Point(5, 1));   Q.push_back(Point(3,  1));
  std::cout << "Q = " << Q << std::endl;

  // Compute the Minkowski sum.
  Polygon_with_holes  sum = CGAL::minkowski_sum_2(P, Q);
  CGAL_assertion(sum.number_of_holes() == 0);
  std::cout << "P (+) Q = " << sum.outer_boundary() << std::endl;
  return 0;
}
