#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPEC;

template <class Kernel>
void test_write_read()
{
  typedef CGAL::Nef_polyhedron_3< Kernel > Nef_polyhedron;
  typedef CGAL::Polyhedron_3< Kernel >     Polyhedron;
  typedef typename Kernel::Point_3         Point;

  typename Kernel::RT n, d;
  std::istringstream str_n("5");
  str_n >> n;
  std::istringstream str_d("1");
  str_d >> d;

  Point p(n, 0, 0, d);
  Point q(0, n, 0, d);
  Point r(0, 0, n, d);
  Point s(0, 0, 0, 1);

  std::cout << "    build...\n";
  Polyhedron P;
  P.make_tetrahedron( p, q, r, s);
  Nef_polyhedron nef_1( P );

  std::cout << "    write...\n";
  std::ofstream out ("temp.nef");
  out << std::setprecision(17) << nef_1;
  out.close();

  std::cout << "    read...\n";
  std::ifstream in ("temp.nef");
  Nef_polyhedron nef_2;
  in >> nef_2;
  in.close();

  std::cout << "    check...\n";
  assert( nef_1 == nef_2);
}

int main()
{
  std::cout << "Testing Exact_predicates_inexact_constructions_kernel\n";
  test_write_read<EPEC>();
  std::cout << "Testing Simple_cartesian<double>\n";
  test_write_read<CGAL::Simple_cartesian<double>>();

  return 0;
}
