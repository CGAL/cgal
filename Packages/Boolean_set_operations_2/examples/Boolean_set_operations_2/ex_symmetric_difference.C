//! \file examples/Boolean_set_operations_2/ex_symmetric_difference.C
// Computing the symmetric difference of two polygons with holes.

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_with_holes_2.h>

#include <list>

typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef Kernel::Point_2                                 Point_2;
typedef CGAL::Polygon_2<Kernel>                         Polygon_2;
typedef CGAL::General_polygon_with_holes_2<Polygon_2>   Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                 Pwh_list_2;

int main ()
{
  // Construct P - a bounded rectangle that contains a rectangular hole.
  Polygon_2   outP;
  Polygon_2   holesP[1];

  outP.push_back (Point_2 (-3, -5));  outP.push_back (Point_2 (3, -5));
  outP.push_back (Point_2 (3, 5));    outP.push_back (Point_2 (-3, 5));
  holesP[0].push_back (Point_2 (-1, -3));
  holesP[0].push_back (Point_2 (-1, 3));
  holesP[0].push_back (Point_2 (1, 3));
  holesP[0].push_back (Point_2 (1, -3));

  Polygon_with_holes_2  P (outP, holesP, holesP + 1); 
  std::cout << "P = " << P << std::endl;

  // Construct Q - a bounded rectangle that contains a rectangular hole.
  Polygon_2   outQ;
  Polygon_2   holesQ[1];

  outQ.push_back (Point_2 (-5, -3));  outQ.push_back (Point_2 (5, -3));
  outQ.push_back (Point_2 (5, 3));    outQ.push_back (Point_2 (-5, 3));
  holesQ[0].push_back (Point_2 (-3, -1));
  holesQ[0].push_back (Point_2 (-3, 1));
  holesQ[0].push_back (Point_2 (3, 1));
  holesQ[0].push_back (Point_2 (3, -1));

  Polygon_with_holes_2  Q (outQ, holesQ, holesQ + 1); 
  std::cout << "Q = " << Q << std::endl;

  // Compute the symmetric difference of P and Q.
  Pwh_list_2                  symmR;
  Pwh_list_2::const_iterator  it;

  CGAL::symmetric_difference (P, Q, std::back_inserter(symmR));

  std::cout << "The symmetric difference:" << std::endl;
  for (it = symmR.begin(); it != symmR.end(); ++it)
    std::cout << "--> " << *it << std::endl;
  
  return (0);
}
