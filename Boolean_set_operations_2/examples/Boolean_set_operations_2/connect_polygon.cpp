/*! \file connect_polygon.cpp
 * Connecting a polygon with holes.
 */

#include "bso_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/connect_holes.h>
#include <list>

typedef CGAL::Cartesian<Number_type>               Kernel;
typedef Kernel::Point_2                            Point_2;
typedef CGAL::Polygon_2<Kernel>                    Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>         Polygon_with_holes_2;

int main (int argc, char **argv)
{
  // Open the input file.
  const char   *filename = "pgn_holes.dat";

  if (argc >= 2)
    filename = argv[1];

  std::ifstream input_file (filename);

  if (! input_file.is_open())
  {
    std::cerr << "Failed to open the " << filename <<std::endl;
    return (1);
  }

  // Read a polygon with holes from a file.
  Polygon_2               outerP;
  unsigned int            num_holes;

  input_file >> outerP;
  input_file >> num_holes;

  std::vector<Polygon_2>  holes (num_holes);
  unsigned int            k;

  for (k = 0; k < num_holes; k++)
    input_file >> holes[k];

  Polygon_with_holes_2    P (outerP, holes.begin(), holes.end());

  // Connect the outer boundary of the polygon with its holes.
  std::list<Point_2>            pts;
  std::list<Point_2>::iterator  pit;

  connect_holes (P, std::back_inserter (pts));

  for (pit = pts.begin(); pit != pts.end(); ++pit)
    std::cout << '(' << *pit << ")  ";
  std::cout << std::endl;

  return (0);
}
