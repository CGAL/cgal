// Test program for Lazy_exact_nt<>.

#include <CGAL/Cartesian.h>
#include <iostream>
#include <CGAL/Random.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Delaunay_triangulation_3.h>

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
typedef leda_real        Exact_NT;
#else
typedef CGAL::Quotient<CGAL::MP_Float> Exact_NT; // doesn't do exact sqrt()
namespace CGAL { Exact_NT sqrt(const Exact_NT &) { abort(); } }
#endif

typedef CGAL::Lazy_exact_nt<Exact_NT> NT;
// typedef Exact_NT NT;

typedef CGAL::Cartesian<NT> Kernel;
typedef Kernel::Point_2     Point;

void predicates()
{
  Point A(NT(1.0)/NT(3),NT(2.0)/NT(3));
  Point B(NT(2.0)/NT(3),NT(3.0)/NT(3));
  Point C(NT(3.0)/NT(3),NT(4.0)/NT(3));
  Point D(NT(4.0)/NT(3),NT(3.0)/NT(3));
  std::cout << "A : " << A << std::endl;
  std::cout << "B : " << B << std::endl;
  std::cout << "C : " << C << std::endl;
  std::cout << (int) CGAL::orientation(A,B,C) << std::endl;
}

namespace CGAL {

template <class NT>
Sign
my_sign(const NT& n)
{
  return CGAL_NTS sign(n);
}

template <class NT>
NT
my_square(const NT& n)
{
  return CGAL_NTS square(n);
}

template <class NT>
NT
my_min(const NT& n, const NT& m)
{
  return min(n, m);
}

} // namespace CGAL

typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >  my_NT;
typedef CGAL::Delaunay_triangulation_3<CGAL::Cartesian<my_NT> > Delaunay;

int my_rand()
{
  return int(CGAL::default_random.get_double()*(1<<31));
}

void delaunay()
{
  Delaunay D;
  for (int i=0; i<100; i++)
    D.insert(Delaunay::Point(my_NT(my_rand()),
                             my_NT(my_rand()),
			     my_NT(my_rand())));
}

int main ()
{
  std::cout.precision(20);
  // NT a;
  // NT b(a);
  // NT c = b;
  NT d (1.0);
  NT e = d + d;
  NT z = CGAL::my_min(e,d);
  (void) d; (void) e; (void) z; // Shut up warnings.
  std::cout << e/NT(3) << std::endl;
  // NT f = abs(NT(0));
  std::cout << "sign(3*(1/3)-1) = " << CGAL::my_sign(NT(3)*(NT(1)/NT(3))-NT(1))
            << std::endl;
  std::cout << "sign(sqrt(2)^2-2) = "
#ifdef CGAL_USE_LEDA
            << CGAL::my_sign(CGAL::my_square(CGAL_NTS sqrt(NT(2)))-NT(2))
#else
            << "EXACT SQUARE ROOT NOT AVAILABLE WITHOUT LEDA"
#endif
            << std::endl;
  // std::cout << "sign(sqrt(2)) = " << CGAL_NTS sign(CGAL_NTS sqrt(NT(2))) << std::endl;
  // std::cout << "sign(square(2)) = " << CGAL_NTS sign(CGAL_NTS square(NT(2))) << std::endl;
  // bool pipo = e < d;
  std::cout << "sizeof(rep) = " << sizeof(CGAL::Lazy_exact_rep<int>)
            << std::endl;
  std::cout << "sizeof(cst) = " << sizeof(CGAL::Lazy_exact_Cst<int>)
            << std::endl;
  std::cout << "sizeof(abs) = " << sizeof(CGAL::Lazy_exact_Abs<int>)
            << std::endl;
  std::cout << "sizeof(add) = " << sizeof(CGAL::Lazy_exact_Add<int>)
            << std::endl;
  std::cout << "sizeof(Exat_NT) = " << sizeof(Exact_NT) << std::endl;
  std::cout << "sizeof(Lazy_exact_nt) = " << sizeof(CGAL::Lazy_exact_nt<int>)
            << std::endl;
  predicates();
  delaunay();
  return 0;
}

