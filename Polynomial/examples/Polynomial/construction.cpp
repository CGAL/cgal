#include <CGAL/config.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  
  typedef CGAL::Polynomial_type_generator<int,2>::Type Poly_2;
  typedef CGAL::Polynomial_traits_d<Poly_2>            PT_2;
  typedef PT_2::Coefficient_type                       Poly_1;
  typedef PT_2::Innermost_coefficient_type             Integer; 
   
  PT_2::Construct_polynomial construct_polynomial;
  
  // constructing a constant polynomial from int
  Poly_2 two(2);
  std::cout << "A constant polynomial: " << two << std::endl;

  
  // construction from an iterator range of univariate polynomials
  
  std::list<Poly_1> univariate_coeffs;
  univariate_coeffs.push_back(Poly_1(3));
  univariate_coeffs.push_back(Poly_1(0));
  univariate_coeffs.push_back(Poly_1(5));
  Poly_2 F = 
    construct_polynomial(univariate_coeffs.begin(),univariate_coeffs.end());
  std::cout << "The bivariate polynomial F: " << F << std::endl;
  
   
  // construction from an iterator range over monomials 
  
  std::list<std::pair<CGAL::Exponent_vector, Integer> > innermost_coeffs;
  CGAL::Exponent_vector ev(std::vector<int>(2,0)); // = [0,0] sequence 
  innermost_coeffs.push_back(std::make_pair(ev,-2));
  ev[0]=3; 
  ev[1]=5; 
  innermost_coeffs.push_back(std::make_pair(ev,2));
  Poly_2 G = 
    construct_polynomial(innermost_coeffs.begin(),innermost_coeffs.end());
  std::cout << "The bivariate polynomial G: " << G << std::endl;
  
  //construction using shift 
  PT_2::Shift shift;
  Poly_2 x = shift(Poly_2(1),1,0); // 'multiply' 1 by x_0^1
  Poly_2 y = shift(Poly_2(1),1,1); // 'multiply' 1 by x_1^1
  
  Poly_2 H = 5 * x * y + 3 * y * y; 
  std::cout << "The bivariate polynomial H: " << H << std::endl;
}
