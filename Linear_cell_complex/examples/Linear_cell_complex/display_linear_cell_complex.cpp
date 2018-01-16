#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_viewer_qt.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;
typedef LCC::Point Point;

int main()
{
  LCC lcc;
  lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                      Point(5,5,0), Point(0,5,0),
                      Point(0,5,4), Point(0,0,4),
                      Point(5,0,4), Point(5,5,4));

  lcc.display_characteristics(std::cout)<<", valid=" 
                                        <<lcc.is_valid()<<std::endl;
  CGAL::draw(lcc);

  return EXIT_SUCCESS;
}
