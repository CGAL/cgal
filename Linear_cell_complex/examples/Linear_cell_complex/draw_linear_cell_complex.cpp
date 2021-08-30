#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;
typedef LCC::Dart_handle Dart_handle;
typedef LCC::Point Point;

int main()
{
  LCC lcc;
  Dart_handle dh1=
    lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                        Point(5,5,0), Point(0,5,0),
                        Point(0,5,4), Point(0,0,4),
                        Point(5,0,4), Point(5,5,4));
  Dart_handle dh2=
    lcc.make_hexahedron(Point(5,0,0), Point(10,0,0),
                        Point(10,5,0), Point(5,5,0),
                        Point(5,5,4), Point(5,0,4),
                        Point(10,0,4), Point(10,5,4));

  lcc.sew<3>(lcc.beta(dh1, 1, 1, 2), lcc.beta(dh2, 2));

  lcc.display_characteristics(std::cout)<<", valid="
                                        <<lcc.is_valid()<<std::endl;
  CGAL::draw(lcc);

  return EXIT_SUCCESS;
}
