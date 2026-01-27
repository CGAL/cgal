#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>

using LCC=CGAL::Linear_cell_complex_for_combinatorial_map<3>;
using Point=LCC::Point;

int main()
{
  LCC lcc;
  lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                      Point(5,5,0), Point(0,5,0),
                      Point(0,5,4), Point(0,0,4),
                      Point(5,0,4), Point(5,5,4));
  CGAL::draw(lcc);

  return EXIT_SUCCESS;
}
