#include <CGAL/Homogeneous.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <fstream>

typedef CGAL::internal::Exact_field_selector<int>::Type   FT;
typedef CGAL::Fraction_traits<FT>::Numerator_type         RT;
typedef CGAL::Lazy_exact_nt < FT >                        Lazy_FT;
typedef CGAL::Lazy_exact_nt < RT >                        Lazy_RT;
typedef CGAL::Homogeneous< Lazy_RT >                      H_Lazy_RT;
typedef CGAL::Simple_homogeneous< Lazy_RT >               SH_Lazy_RT;


template <class Kernel>
void test_write_read()
{
  typedef CGAL::Nef_polyhedron_3< Kernel > Nef_polyhedron;
  typedef CGAL::Polyhedron_3< Kernel >     Polyhedron;
  typedef typename Kernel::Point_3         Point;

  typename Kernel::RT n( std::string("6369051672525773"));
  typename Kernel::RT d( std::string("4503599627370496"));

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
  std::cout << "Testing Homogeneous< Lazy_exact_nt<RT> >\n";
  test_write_read<H_Lazy_RT>();
  std::cout << "Testing Simple_homogeneous< Lazy_exact_nt<RT> >\n";
  test_write_read<SH_Lazy_RT>();

  return 0;
}
