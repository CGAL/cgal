/*! \file do_intersect.cpp
 * Determining whether two triangles intersect.
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

#include "print_utils.h"

int main() {
  Polygon_2 P;
  P.push_back(Point_2(-1, 1));
  P.push_back(Point_2(0, -1));
  P.push_back(Point_2(1, 1));
  std::cout << "P = "; print_polygon(P);

  Polygon_2 Q;
  Q.push_back(Point_2(-1, -1));
  Q.push_back(Point_2(1, -1));
  Q.push_back(Point_2(0, 1));
  std::cout << "Q = "; print_polygon(Q);

  if ((CGAL::do_intersect(P, Q)))
    std::cout << "The two polygons intersect in their interior." << std::endl;
  else
    std::cout << "The two polygons do not intersect." << std::endl;

  return 0;
}
