#include <CGAL/Complex_without_sqrt.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Interval_nt.h>
#include <iostream>
#include <sstream>

using namespace CGAL;

typedef Complex_without_sqrt<Gmpq> Complex_gmpq;
typedef Complex_without_sqrt<Exact_integer> Complex_integer;
typedef Complex_without_sqrt<Interval_nt<>> Complex_interval;


int main() {
  // Complex_gmpq tests :
  Complex_gmpq zero_gmpq = Complex_gmpq ();
  assert( zero_gmpq == Complex_gmpq(Gmpq(0), Gmpq(0)) );

  Complex_gmpq one_gmpq (Gmpq(1));
  assert( one_gmpq == Complex_gmpq(Gmpq(1), Gmpq(0)) );

  Complex_gmpq z1_gmpq (Gmpq(1,2), Gmpq(-3));
  z1_gmpq = -z1_gmpq;

  Complex_gmpq z2_gmpq;
  z2_gmpq.set_real_part(Gmpq(-5,7));
  z2_gmpq.set_imaginary_part(Gmpq(11,13));
  z2_gmpq = z2_gmpq.conjugate() + z1_gmpq - one_gmpq;

  assert( - z1_gmpq * z1_gmpq / z2_gmpq == -Complex_gmpq(Gmpq(855491,632146), Gmpq(844298,316073)) );

  assert( z1_gmpq.real_part() == Gmpq(-1,2) );
  assert( z1_gmpq.imaginary_part() == Gmpq(3) );
  assert( z1_gmpq.squared_modulus() == Gmpq(37,4) );
  assert( z1_gmpq != z2_gmpq);
  assert( z2_gmpq == z2_gmpq );
  assert( z2_gmpq == Complex_gmpq(Gmpq(-31,14), Gmpq(28,13)) );

  std::cout << "printing a complex for test purposes : " << std::endl << z2_gmpq << std::endl;

  Complex_gmpq z3_gmpq;
  std::stringstream buffer;
  buffer << z2_gmpq;
  buffer >> z3_gmpq;
  assert( z3_gmpq == z2_gmpq );

  // Complex_integer tests :
  Complex_integer zero_integer = Complex_integer ();
  assert( zero_integer == Complex_integer(Exact_integer(0),Exact_integer(0)) );

  Complex_integer one_integer (Exact_integer(1));
  assert( one_integer == Complex_integer(Exact_integer(1), Exact_integer(0)) );

  Complex_integer z1_integer (Exact_integer(17), Exact_integer(-13));
  z1_integer = -z1_integer;

  Complex_integer z2_integer;
  z2_integer.set_real_part(Exact_integer(-7));
  z2_integer.set_imaginary_part(Exact_integer(43));
  z2_integer = z2_integer.conjugate() + z1_integer - one_integer;

  assert( z1_integer * z1_integer / z2_integer == Complex_integer(Exact_integer(0),Exact_integer(0)) );

  assert( z1_integer.real_part() ==  Exact_integer(-17) );
  assert( z1_integer.imaginary_part() ==  Exact_integer(13) );
  assert( z1_integer.squared_modulus() == 458 );
  assert( z1_integer != z2_integer);
  assert( z2_integer == z2_integer );
  assert( z2_integer == Complex_integer(Exact_integer(-25),Exact_integer(-30)) );

  std::cout << "printing a complex for test purposes : " << std::endl << z2_integer << std::endl;

  Complex_integer z3_integer;
  buffer << z2_integer;
  buffer >> z3_integer;
  assert( z3_integer == z2_integer );

  // Complex_interval test :
  Complex_interval z_interval (Interval_nt<>(1, 2), Interval_nt<>(1, 2));
  assert( (z_interval * z_interval / Complex_interval(Interval_nt<>(5, 6))).squared_modulus() < Interval_nt<>(10,20) );

  return 0;
}
