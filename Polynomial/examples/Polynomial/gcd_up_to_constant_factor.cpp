#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  typedef CGAL::Polynomial_type_generator<int,1>::Type Poly_1;
  typedef CGAL::Polynomial_traits_d<Poly_1>            PT_1;

  PT_1::Shift                     shift;
  PT_1::Gcd                       gcd;
  PT_1::Gcd_up_to_constant_factor gcd_utcf;
  PT_1::Multivariate_content      mcontent;
  PT_1::Canonicalize              canonicalize;

  //construction using shift
  Poly_1 x = shift(Poly_1(1),1,0); // x^1

  // common factor 7 * (x^2-2)
  Poly_1 F = 21*(x-5)*(x*x-2); // = 21*x^3 + (-105)*x^2 + (-42)*x + 210
  Poly_1 G = 14*(x-3)*(x*x-2); // = 14*x^3 + (-42)*x^2 + (-28)*x + 84

  std::cout << "The univariate polynomial F: " << F << std::endl;
  std::cout << "The univariate polynomial G: " << G << std::endl;
  std::cout << "Common multivariate content:              "
            << CGAL::gcd(mcontent(F),mcontent(G)) // = 7
            << std::endl;
  std::cout << "The gcd of F and G:                       "
            << gcd(F,G)                           // = 7*x^2 + (-14)
            << std::endl;
  std::cout << "The gcd up to constant factor of F and G: "
            << gcd_utcf(F,G)                      // = x^2 + (-2)
            << std::endl;
  std::cout << "Same as canonicalized gcd of F and G:     "
            << canonicalize(gcd_utcf(F,G))        // = x^2 + (-2)
            << std::endl;

}
