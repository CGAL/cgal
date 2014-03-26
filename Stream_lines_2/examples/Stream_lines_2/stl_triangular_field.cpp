#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Stream_lines_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Triangular_field_2.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Point_2                                                  Point;
typedef K::Vector_2                                                 Vector;
typedef CGAL::Triangular_field_2<K>                                 Field;
typedef CGAL::Runge_kutta_integrator_2<Field>                       Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>         Strl;
typedef Strl::Stream_line_iterator_2                                stl_iterator;

int main()
{
  Runge_kutta_integrator runge_kutta_integrator(1);

  /*datap.tri.cin and datav.tri.cin are ascii files where are stored the vector values*/
  std::ifstream inp("data/datap.tri.cin");
  std::ifstream inv("data/datav.tri.cin");
  std::istream_iterator<Point> beginp(inp);
  std::istream_iterator<Vector> beginv(inv);
  std::istream_iterator<Point> endp;

  Field triangular_field(beginp, endp, beginv);

  /* the placement of streamlines */
  std::cout << "processing...\n";
  double dSep = 30.0;
  double dRat = 1.6;
  Strl Stream_lines(triangular_field, runge_kutta_integrator,dSep,dRat);
  std::cout << "placement generated\n";

  /*writing streamlines to streamlines.stl */
  std::cout << "streamlines.stl\n";
  std::ofstream fw("streamlines.stl",std::ios::out);
  Stream_lines.print_stream_lines(fw);
}
