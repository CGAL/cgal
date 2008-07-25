#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Gmpz.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  
  typedef CGAL::Polynomial_type_generator<CGAL::Gmpz,2>::Type Poly_int_2;

  typedef CGAL::Polynomial_traits_d<Poly_int_2>            PT_2;
  typedef PT_2::Coefficient                                Poly_int_1;
  typedef PT_2::Innermost_coefficient                      Integer; // this is Gmpz 
  
  // constructing a constant polynomial from int
  Poly_int_2 two(2);
  std::cout << "A constant polynomial: " << two << std::endl;

  // construction from an iterator range of univariate polynomials
  std::list<Poly_int_1> univariate_coeffs;
  univariate_coeffs.push_back(Poly_int_1(3));
  univariate_coeffs.push_back(Poly_int_1(0));
  univariate_coeffs.push_back(Poly_int_1(5));
  
  Poly_int_2 F = 
  PT_2::Construct_polynomial()(univariate_coeffs.begin(),univariate_coeffs.end());
  std::cout << "The bivariate polynomial F: " << F << std::endl;
  
  // construction from an iterator range over monomials 
  std::list<std::pair<CGAL::Exponent_vector, Integer> > innermost_coeffs;
  CGAL::Exponent_vector ev(2,0); // = [0,0] sequenze 
  innermost_coeffs.push_back(std::make_pair(ev,CGAL::Gmpz(-2)));
  ev[0]=3; ev[1]=5; 
  innermost_coeffs.push_back(std::make_pair(ev,CGAL::Gmpz(2)));
  Poly_int_2 G = 
    PT_2::Construct_polynomial()(innermost_coeffs.begin(),innermost_coeffs.end());
  std::cout << "The bivariate polynomial G: " << G << std::endl;
  
  //construction using shift 
  Poly_int_2 x = PT_2::Shift()(Poly_int_2(1),1,0); // shift 1 by x_0^1
  Poly_int_2 y = PT_2::Shift()(Poly_int_2(1),1,1); // shift 1 by x_1^1
  
  Poly_int_2 H = 5 * x * y + 3 * CGAL::ipower(y,5); 
  std::cout << "The bivariate polynomial H: " << H << std::endl;
}
