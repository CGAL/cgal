// $URL$
// $Id$

#include <CGAL/config.h>
#ifdef CGAL_USE_MPFI
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Gmpz.h>
#include <vector>

typedef CGAL::Algebraic_kernel_d_1<CGAL::Gmpz>          AK;
typedef AK::Coefficient                                 Coefficient;
typedef AK::Polynomial_1                                Polynomial_1;
typedef AK::Algebraic_real_1                            Algebraic_real_1;
typedef AK::Bound                                       Bound;
typedef std::pair<Bound,Bound>                          Interval;

int main(){
  AK ak;

  AK::Construct_algebraic_real_1 construct_algebraic_real_1 = ak.construct_algebraic_real_1_object();
  Polynomial_1 x = CGAL::shift(AK::Polynomial_1(1),1); // the monomial x
  Algebraic_real_1 a = construct_algebraic_real_1(x*x-2,1); //  sqrt(2)
  Algebraic_real_1 b = construct_algebraic_real_1(x*x-3,1); //  sqrt(3)

  // Algebraic_real_1 is RealEmbeddable (just some functions:)
  std::cout << "sign of a is                 : " << CGAL::sign(a)      << "\n";
  std::cout << "double approximation of a is : " << CGAL::to_double(a) << "\n";
  std::cout << "double approximation of b is : " << CGAL::to_double(b) << "\n";
  std::cout << "double lower bound of a      : " << CGAL::to_interval(a).first  << "\n";
  std::cout << "double upper bound of a      : " << CGAL::to_interval(a).second << "\n";
  std::cout << "LessThanComparable (a<b)     : " << (a<b) << "\n\n";

  // use compare_1 with int, Bound, Coefficient, Algebraic_real_1
  AK::Compare_1 compare_1 = ak.compare_1_object();
  std::cout << " compare with an int                  : " << compare_1(a ,int(2)) << "\n";
  std::cout << " compare with an Coefficient          : " << compare_1(a ,Coefficient(2)) << "\n";
  std::cout << " compare with an Bound                : " << compare_1(a ,Bound(2)) << "\n";
  std::cout << " compare with another Algebraic_real_1: " << compare_1(a ,b) << "\n\n";

  // get a value between two roots
  AK::Bound_between_1 bound_between_1 = ak.bound_between_1_object();
  std::cout << " value between sqrt(2) and sqrt(3) " << bound_between_1(a,b) << "\n";
  std::cout << " is larger than sqrt(2)            " << compare_1(bound_between_1(a,b),a) << "\n";
  std::cout << " is less   than sqrt(3)            " << compare_1(bound_between_1(a,b),b) << "\n\n";

  // approximate with relative precision
  AK::Approximate_relative_1 approx_r = ak.approximate_relative_1_object();
  std::cout << " lower bound of a with at least 100 bits:    "<< approx_r(a,100).first  << "\n";
  std::cout << " upper bound of a with at least 100 bits:    "<< approx_r(a,100).second << "\n\n";

  // approximate with absolute error
  AK::Approximate_absolute_1 approx_a = ak.approximate_absolute_1_object();
  std::cout << " lower bound of b with error less than 2^-100:   "<< approx_a(b,100).first  << "\n";
  std::cout << " upper bound of b with error less than 2^-100:   "<< approx_a(b,100).second << "\n\n";

  return 0;
  }
#else
int main(){
  std::cout << "This example requires CGAL to be configured with library MPFI." << std::endl;
return 0;
}
#endif
