/*! \file sequence.cpp
 * Performing a sequence of Boolean set-operations.
 */

#include "bso_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>

#include <list>

// instead of
//typedef CGAL::Cartesian<Number_type>               Kernel;
// workaround for VC++
struct Kernel : public CGAL::Cartesian<Number_type> {};

typedef Kernel::Point_2                            Point_2;
typedef CGAL::Polygon_2<Kernel>                    Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>         Polygon_with_holes_2;
typedef CGAL::Polygon_set_2<Kernel>                Polygon_set_2;

#include "print_utils.h"

int main ()
{
  // Construct the two initial polygons and the clipping rectangle.
  Polygon_2 P;
  P.push_back (Point_2 (0, 1));
  P.push_back (Point_2 (2, 0));
  P.push_back (Point_2 (1, 1));
  P.push_back (Point_2 (2, 2));

  Polygon_2 Q;
  Q.push_back (Point_2 (3, 1));
  Q.push_back (Point_2 (1, 2));
  Q.push_back (Point_2 (2, 1));
  Q.push_back (Point_2 (1, 0));

  Polygon_2 rect;
  rect.push_back (Point_2 (0, 0));
  rect.push_back (Point_2 (3, 0));
  rect.push_back (Point_2 (3, 2));
  rect.push_back (Point_2 (0, 2));

  // Perform a sequence of operations.
  Polygon_set_2 S;
  S.insert (P);
  S.join (Q);                   // Compute the union of P and Q.
  S.complement();               // Compute the complement.
  S.intersection (rect);        // Intersect with the clipping rectangle.

  // Print the result.
  std::list<Polygon_with_holes_2> res;
  std::list<Polygon_with_holes_2>::const_iterator it;

  std::cout << "The result contains " << S.number_of_polygons_with_holes()
            << " components:" << std::endl;

  S.polygons_with_holes (std::back_inserter (res));
  for (it = res.begin(); it != res.end(); ++it) {
    std::cout << "--> ";
    print_polygon_with_holes (*it);
  }

  return 0;
}
