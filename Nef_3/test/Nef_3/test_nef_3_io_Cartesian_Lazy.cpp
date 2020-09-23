#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <fstream>

typedef CGAL::internal::Exact_field_selector<int>::Type   FT;
typedef CGAL::Lazy_exact_nt < FT >                        Lazy_FT;
typedef CGAL::Cartesian< Lazy_FT >                        C_Lazy_FT;
typedef CGAL::Simple_cartesian< Lazy_FT >                 SC_Lazy_FT;


template <class Kernel>
void test_write_read()
{
  typedef CGAL::Nef_polyhedron_3< Kernel > Nef_polyhedron;
  typedef CGAL::Polyhedron_3< Kernel >     Polyhedron;
  typedef typename Kernel::Point_3         Point;
  typedef typename Kernel::RT              RT;
  // Boost.Multiprecision has trouble constructing Lazy_FT directly from string
  typedef typename RT::Exact_type          RTE;

  RT n(RTE(std::string("6369051672525773")));
  RT d(RTE(std::string("4503599627370496")));

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
  out << nef_1;
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
  std::cout << "Testing Cartesian< Lazy_exact_nt<FT> >\n";
  test_write_read<C_Lazy_FT>();
  std::cout << "Testing Simple_cartesian< Lazy_exact_nt<FT> >\n";
  test_write_read<SC_Lazy_FT>();

  return 0;
}
