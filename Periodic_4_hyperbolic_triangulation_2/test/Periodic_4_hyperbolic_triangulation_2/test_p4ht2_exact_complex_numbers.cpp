#include <CGAL/Exact_algebraic.h>
#include <CGAL/internal/Exact_complex.h>
#include <CGAL/Cartesian.h>

#include <iostream>

int main(int, char**)
{
  typedef CGAL::Exact_algebraic       NT;
  typedef CGAL::Cartesian<NT>         Kernel;
  typedef CGAL::Exact_complex<NT>     cplx;
  typedef Kernel::Point_2             Point;
  typedef Kernel::Circle_2            Circle;

  cplx n;
  n.set_imag(-1);
  n.set_real(3);
  std::cout << "n = " << n << std::endl;
  std::cout << "sq_mod(n) = " << n.square_modulus() << std::endl;
  std::cout << "mod(n)    = " << n.modulus() << std::endl;
  std::cout << std::endl;
  cplx m(2, 3);
  std::cout << "m = " << m << std::endl;
  std::cout << "sq_mod(m) = " << m.square_modulus() << std::endl;
  std::cout << "mod(m)    = " << m.modulus() << std::endl;
  std::cout << std::endl;

  std::cout << "m + n = " << m + n << std::endl;
  std::cout << "m - n = " << m - n << std::endl;
  std::cout << "m * n = " << m * n << std::endl;
  std::cout << "m / n = " << m / n << std::endl;
  std::cout << "n / m = " << n / m << std::endl;

  cplx mon = m/n;
  cplx nom = n/m;

  cplx res1(NT(3)/NT(10), NT(11)/NT(10));
  cplx res2(NT(3)/NT(13), NT(-11)/NT(13));

  assert(mon == res1);
  assert(nom == res2);

  std::cout << "n < m: " << (n < m ? "true" : "false") << std::endl;

  Point pt(NT(2)/NT(10), NT(6)/NT(10));
  cplx p(pt.x(), pt.y());
  std::cout << "p = " << p << std::endl;
  std::cout << "recip = " << p.reciprocal() << std::endl;
  std::cout << "inver = " << p.invert_in_unit_circle() << std::endl;

  Circle c(Point(NT(2)/NT(10), NT(3)/NT(4)), NT(1)/NT(2));
  cplx ipc = p.invert_in_circle(c);
  std::cout << "inverted in circle: " << ipc << std::endl;

  return EXIT_SUCCESS;
}
