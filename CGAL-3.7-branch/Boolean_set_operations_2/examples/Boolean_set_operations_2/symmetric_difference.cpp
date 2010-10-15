/*! \file symmetric_difference.cpp
 * Computing the symmetric difference of two polygons with holes.
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <list>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                   Pwh_list_2;

#include "print_utils.h"

int main ()
{
  // Construct P - a bounded rectangle that contains a rectangular hole.
  Polygon_2 outP;
  Polygon_2 holesP[1];

  outP.push_back (Point_2 (-3, -5));  outP.push_back (Point_2 (3, -5));
  outP.push_back (Point_2 (3, 5));    outP.push_back (Point_2 (-3, 5));
  holesP[0].push_back (Point_2 (-1, -3));
  holesP[0].push_back (Point_2 (-1, 3));
  holesP[0].push_back (Point_2 (1, 3));
  holesP[0].push_back (Point_2 (1, -3));

  Polygon_with_holes_2  P (outP, holesP, holesP + 1);
  std::cout << "P = "; print_polygon_with_holes (P);

  // Construct Q - a bounded rectangle that contains a rectangular hole.
  Polygon_2 outQ;
  Polygon_2 holesQ[1];

  outQ.push_back (Point_2 (-5, -3));  outQ.push_back (Point_2 (5, -3));
  outQ.push_back (Point_2 (5, 3));    outQ.push_back (Point_2 (-5, 3));
  holesQ[0].push_back (Point_2 (-3, -1));
  holesQ[0].push_back (Point_2 (-3, 1));
  holesQ[0].push_back (Point_2 (3, 1));
  holesQ[0].push_back (Point_2 (3, -1));

  Polygon_with_holes_2  Q (outQ, holesQ, holesQ + 1);
  std::cout << "Q = "; print_polygon_with_holes (Q);

  // Compute the symmetric difference of P and Q.
  Pwh_list_2 symmR;
  Pwh_list_2::const_iterator it;

  CGAL::symmetric_difference (P, Q, std::back_inserter(symmR));

  std::cout << "The symmetric difference:" << std::endl;
  for (it = symmR.begin(); it != symmR.end(); ++it) {
    std::cout << "--> ";
    print_polygon_with_holes (*it);
  }

  return 0;
}
