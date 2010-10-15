// $URL$
// $Id$

#include <CGAL/basic.h>
#ifdef CGAL_USE_MPFI 
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Gmpz.h>
#include <vector>

typedef CGAL::Algebraic_kernel_d_1<CGAL::Gmpz>          AK;
typedef AK::Polynomial_1                                Polynomial_1;
typedef AK::Algebraic_real_1                            Algebraic_real_1;
typedef AK::Coefficient                                 Coefficient;
typedef AK::Bound                                       Bound;
typedef AK::Multiplicity_type                           Multiplicity_type;

int main(){
  AK ak; // an object of 
  AK::Construct_algebraic_real_1 construct_algreal_1 = ak.construct_algebraic_real_1_object();
  AK::Isolate_1 isolate_1 = ak.isolate_1_object();
  AK::Compute_polynomial_1 compute_polynomial_1 = ak.compute_polynomial_1_object();

  // construct an algebraic number from an integer
  Algebraic_real_1 frominteger=construct_algreal_1(int(2));
  std::cout << "Construct from int: " << frominteger << "\n";

  // the constructed algebraic number is root of a polynomial
  Polynomial_1 pol=compute_polynomial_1(frominteger);
  std::cout << "The constructed number is root of: " << pol << "\n";

  // construct an algebraic number from a polynomial and an isolating interval
  Polynomial_1 x = CGAL::shift(AK::Polynomial_1(1),1); // the monomial x
  Algebraic_real_1 frominterval=construct_algreal_1(x*x-2,Bound(0),Bound(2));
  std::cout << "Construct from isolating interval: " << frominterval << "\n";

  // isolate the second algebraic number from the first: this is to say,
  // isolating the second algebraic number with respect to the polynomial
  // of which the first constructed number is root
  std::pair<Bound,Bound> isolation1 = isolate_1(frominterval,pol);
  std::cout << "Isolating the second algebraic number gives: ["
            << isolation1.first << "," << isolation1.second << "]\n";

  // isolate again the same algebraic number, this time with respect to
  // the polynomial 10*x-14 (which has root 1.4, close to this algebraic
  // number)
  std::pair<Bound,Bound> isolation2 = isolate_1(frominterval,10*x-14);
  std::cout << "Isolating again the second algebraic number gives: ["
            << isolation2.first << "," << isolation2.second << "]\n";

  return 0;
}
#else
int main(){
  std::cout << "This example requires CGAL to be configured with library MPFI." << std::endl;
return 0;
}
#endif
