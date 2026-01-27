#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/draw_linear_cell_complex.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC1;
typedef CGAL::Linear_cell_complex_for_generalized_map<3> LCC2;

template<typename LCC>
void test()
{
  LCC lcc;
  using Point=typename LCC::Point;

  typename LCC::Dart_descriptor d1=
    lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                        Point(5,5,0), Point(0,5,0),
                        Point(0,5,4), Point(0,0,4),
                        Point(5,0,4), Point(5,5,4));
  typename LCC::Dart_descriptor d2=
    lcc.make_quadrangle(Point(5,2,2), Point(5,1,2),
                        Point(5,1,1), Point(5,2,1));

  lcc.insert_cell_1_between_two_cells_2
    (lcc.template opposite<2>(lcc.next(lcc.next(d1))),
     lcc.next(lcc.next(d2)));

  CGAL::draw(lcc);
}

int main()
{
  test<LCC1>();
  test<LCC2>();
  return EXIT_SUCCESS;
}
