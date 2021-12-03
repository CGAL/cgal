/*! \file connect_polygon.cpp
 * Connecting a polygon with holes.
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/connect_holes.h>
#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;

int main (int argc, char* argv[])
{

  // Get the name of the input file from the command line, or use the default
  // pgn_holes.dat file if no command-line parameters are given.
  //more data files can be found under test data
  //boundary no other connections are made.
  const char* filename = (argc > 1) ? argv[1] : "pgn_holes.dat";
  std::ifstream input_file (filename);
  if (! input_file.is_open())
  {
    std::cerr << "Failed to open the " << filename <<std::endl;
    return -1;
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
