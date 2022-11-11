// $URL$
// $Id$

#include <CGAL/config.h>
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
  AK ak;
  AK::Solve_1 solve_1 = ak.solve_1_object();
  AK::Sign_at_1 sign_at_1 = ak.sign_at_1_object();
  AK::Is_zero_at_1 is_zero_at_1 = ak.is_zero_at_1_object();

  // construct the polynomials p=x^2-5 and q=x-2
  Polynomial_1 x = CGAL::shift(AK::Polynomial_1(1),1); // the monomial x
  Polynomial_1 p = x*x-5;
  std::cout << "Polynomial p: " << p << "\n";
  Polynomial_1 q = x-2;
  std::cout << "Polynomial q: " << q << "\n";

  // find the roots of p (it has two roots) and q (one root)
  std::vector<Algebraic_real_1> roots_p,roots_q;
  solve_1(p,true, std::back_inserter(roots_p));
  solve_1(q,true, std::back_inserter(roots_q));

  // evaluate the second root of p in q
  std::cout << "Sign of the evaluation of root 2 of p in q: "
            << sign_at_1(q,roots_p[1]) << "\n";

  // evaluate the root of q in p
  std::cout << "Sign of the evaluation of root 1 of q in p: "
            << sign_at_1(p,roots_q[0]) << "\n";

  // check whether the evaluation of the first root of p in p is zero
  std::cout << "Is zero the evaluation of root 1 of p in p? "
            << is_zero_at_1(p,roots_p[0]) << "\n";

  return 0;
}

#else
int main(){
  std::cout << "This example requires CGAL to be configured with library MPFI." << std::endl;
return 0;
}
#endif
