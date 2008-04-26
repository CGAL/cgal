//! \file examples/Minkowski_sum_2/sum_triangles.cpp
// Computing the Minkowski sum of two triangles.

#include "ms_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/minkowski_sum_2.h>
#include <iostream>

#include "print_utils.h"

typedef CGAL::Cartesian<Number_type>                Kernel;
typedef Kernel::Point_2                             Point_2;
typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

int main ()
{
  // Construct the first polygon (a triangle).
  Polygon_2   P;

  P.push_back (Point_2 (0, 0));
  P.push_back (Point_2 (6, 0));
  P.push_back (Point_2 (3, 5));

  // Construct the second polygon (a triangle).
  Polygon_2   Q;

  Q.push_back (Point_2 (0, 0));
  Q.push_back (Point_2 (2, -2));
  Q.push_back (Point_2 (2, 2));

  // Compute the Minkowski sum.
  Polygon_with_holes_2  sum = minkowski_sum_2 (P, Q);

  CGAL_assertion (sum.number_of_holes() == 0);

  std::cout << "P = "; print_polygon (P);
  std::cout << "Q = "; print_polygon (Q);
  std::cout << "P (+) Q = "; print_polygon (sum.outer_boundary());

  return (0);
}
