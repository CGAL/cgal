// Test program for Lazy_exact_nt<>.

#include <CGAL/Cartesian.h>
#include <CGAL/use.h>
#include <iostream>
#include <cassert>
#include <cstdlib>

#include <CGAL/Random.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Delaunay_triangulation_3.h>


#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
typedef leda_real        Exact_NT;
//#elif defined CGAL_USE_GMP
//#  include <CGAL/Gmpq.h>
//typedef CGAL::Gmpq       Exact_NT;
#else
typedef CGAL::Quotient<CGAL::MP_Float> Exact_NT; // doesn't do exact sqrt()
namespace CGAL {
Exact_NT sqrt(const Exact_NT &)
{
  std::abort();
  return Exact_NT();
}
}
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
#ifdef CGAL_PROFILE
  assert(A.x().depth() == 1);
#endif
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
    return ::CGAL::min BOOST_PREVENT_MACRO_SUBSTITUTION (n, m);
}

} // namespace CGAL

typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >  my_NT;
typedef CGAL::Delaunay_triangulation_3<CGAL::Cartesian<my_NT> > Delaunay;

int my_rand()
{
  return int(CGAL::get_default_random().get_double()*(1<<31));
}

void delaunay()
{
  Delaunay D;
  for (int i=0; i<100; i++)
    D.insert(Delaunay::Point(my_NT(my_rand()),
                             my_NT(my_rand()),
			     my_NT(my_rand())));
}

// Tests the precision of to_double()
void test_to_double()
{
  std::cout << "Current relative precision of to_double is : "
            << NT::get_relative_precision_of_to_double() << std::endl;
  double prec = 1.0/(1<<30)/(1<<10);
  std::cout << "Setting it to : " << prec << std::endl;
  NT::set_relative_precision_of_to_double(prec);
  assert(NT::get_relative_precision_of_to_double() == prec);

  // First compute an approximated value for 1.
  NT one = 1;
  NT three = 3;
  NT tmp = one/three;
  tmp *= three;
  tmp = (tmp + 1/tmp)/2;
  std::cout << "Approximated interval for 1 : " << tmp.approx() << std::endl;

  // Now we square it repeatedly (the interval is going to grow), and we check
  // that to_double() stays reasonnably close to 1.
  for (int i = 0; i < 20; ++i) {
    tmp = CGAL_NTS square(tmp);
    double d = CGAL_NTS to_double(tmp);
    std::cout << "double approximation is : " << d << std::endl;
    std::cout << "interval approximation is : " << tmp.approx() << std::endl;
    //std::cout << "numerator   = " << tmp.exact().numerator() << std::endl;
    //std::cout << "denominator = " << tmp.exact().denominator() << std::endl;
    assert(d > 0.999 && d < 1.001);
  }
}

int main ()
{
  std::cout.precision(20);
  // NT a;
  // NT b(a);
  // NT c = b;
  NT d (1.0);
  NT e = d + d;
  NT f = (short)1;
  NT g = CGAL::POSITIVE;
  NT z = CGAL::my_min(e,d);
  (void) d; (void) e; (void) f; (void) g; (void) z; // Shut up warnings.
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
  std::cout << "sizeof(Exat_NT) = " << sizeof(Exact_NT) << std::endl;
  std::cout << "sizeof(Lazy_exact_nt) = " << sizeof(CGAL::Lazy_exact_nt<int>)
            << std::endl;
  predicates();
  delaunay();

  // Try conversions from Lazy_exact_nt<XXX> of different XXXs.
  CGAL::Lazy_exact_nt<CGAL::MP_Float> one(1), two;
  two = one + one;
  CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > eins(1), zwei;
  zwei = eins + eins;
  CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > deux(two);
  assert(zwei == deux);
#ifdef CGAL_PROFILE
  assert(zwei.depth() == 1);
#endif

  test_to_double();

  // Test % and gcd.
  CGAL::Lazy_exact_nt<int> tmp35(35), tmp7(7), tmp9(9);
  assert((tmp35 % tmp7) == 0);
  assert((tmp35 % 7) == 0);
  assert(CGAL_NTS gcd(tmp35, tmp7) == 7);
  tmp9 %= tmp7;
  assert(tmp9 == 2);
  tmp9 = 9;
  tmp9 %= 7;
  assert(tmp9 == 2);

  // Test if all operators are defined.
  int r = 1;
  CGAL::Lazy_exact_nt<double> q(1);

  q+q; q+r; r+q; q+1; 1+q;
  q-q; q-r; r-q; q-1; 1-q;
  q*q; q*r; r*q; q*1; 1*q;
  q/q; q/r; r/q; q/1; 1/q;
  -q;
  bool b; // avoid clang warning: equality comparison result unused [-Wunused-comparison]
  b = q<q; b = q<r; b = r<q; b = q<1; b = 1<q;
  b = q>q; b = q>r; b = r>q; b = q>1; b = 1>q;
  b = q<=q; b = q<=r; b = r<=q; b = q<=1; b = 1<=q;
  b = q>=q; b = q>=r; b = r>=q; b = q>=1; b = 1>=q;

  b = q==q; b = q==r; b = r==q; b = q==1; b = 1==q;
  b = q!=q; b = q!=r; b = r!=q; b = q!=1; b = 1!=q;

  // Test comparisons with double.
  b = two<1.0; b = 1.0<two;
  b = two>1.0; b = 1.0>two;
  b = two<=1.0; b = 1.0<=two;
  b = two>=1.0; b = 1.0>=two;
  b = two!=1.0; b = 1.0!=two;
  b = two==1.0; b = 1.0==two;

  b = zwei<1.0; b = 1.0<zwei;
  b = zwei>1.0; b = 1.0>zwei;
  b = zwei<=1.0; b = 1.0<=zwei;
  b = zwei>=1.0; b = 1.0>=zwei;
  b = zwei!=1.0; b = 1.0!=zwei;
  b = zwei==1.0; b = 1.0==zwei;

  CGAL_USE(b);

  return 0;
}
