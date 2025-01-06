#include <CGAL/Complex_number.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>
#include <iostream>
#include <sstream>

typedef CGAL::Exact_rational Exact_rational;
typedef CGAL::Interval_nt<> Interval;
typedef CGAL::Complex_number<Exact_rational> Complex_rational;
typedef CGAL::Complex_number<Interval> Complex_interval;

int main() {
  // Complex_rational tests :
  Complex_rational zero_rational = Complex_rational ();
  assert( zero_rational == Complex_rational(Exact_rational(0), Exact_rational(0)) );

  Complex_rational one_rational (Exact_rational(1));
  assert( one_rational == Complex_rational(Exact_rational(1), Exact_rational(0)) );

  Complex_rational z1_rational (Exact_rational(1,2), Exact_rational(-3));
  z1_rational = - z1_rational;

  Complex_rational z2_rational;
  z2_rational.real(Exact_rational(-5,7));
  z2_rational.imag(Exact_rational(11,13));
  z2_rational = conj(z2_rational) + z1_rational - one_rational;

  assert( - z1_rational * z1_rational / z2_rational == -Complex_rational(Exact_rational(855491,632146), Exact_rational(844298,316073)) );
  assert( z1_rational.real() == Exact_rational(-1,2) );
  assert( z1_rational.imag() == Exact_rational(3) );
  assert( norm(z1_rational) == Exact_rational(37,4) );
  assert( z1_rational != z2_rational);
  assert( z2_rational == z2_rational );
  assert( z2_rational == Complex_rational(Exact_rational(-31,14), Exact_rational(28,13)) );

  std::cout << "printing a complex for test purposes : " << std::endl << z2_rational << std::endl;

  Complex_rational z3_rational;
  std::stringstream buffer;
  buffer << z2_rational;
  buffer >> z3_rational;
  assert( z3_rational == z2_rational );

  // Complex_interval test :
  Complex_interval z_interval (Interval(1, 2), Interval(1, 2));
  assert( norm(z_interval * z_interval / Complex_interval(Interval(5, 6))) < Interval(10,20) );

  return 0;
}
