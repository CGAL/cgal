#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  typedef CGAL::Polynomial_type_generator<int,2>::Type Poly_2;
  typedef CGAL::Polynomial_traits_d<Poly_2>            PT_2;
  
  //construction using shift 
  Poly_2 x = PT_2::Shift()(Poly_2(1),1,0); // = x^1
  Poly_2 y = PT_2::Shift()(Poly_2(1),1,1); // = y^1
  
  Poly_2 F // = (11*x^2 + 5*x)*y^4 + (7*x^2)*y^3
    = 11 * CGAL::ipower(y,4) * CGAL::ipower(x,2) 
    + 5 * CGAL::ipower(y,4)  * CGAL::ipower(x,1) 
    + 7 * CGAL::ipower(y,3)  * CGAL::ipower(x,2);  
  std::cout << "The bivariate polynomial F: " << F <<"\n"<< std::endl;
  
  PT_2::Get_coefficient get_coefficient;
  std::cout << "Coefficient of y^0: "<< get_coefficient(F,0) << std::endl;
  std::cout << "Coefficient of y^1: "<< get_coefficient(F,1) << std::endl;
  std::cout << "Coefficient of y^2: "<< get_coefficient(F,2) << std::endl;
  std::cout << "Coefficient of y^3: "<< get_coefficient(F,3) << std::endl;
  std::cout << "Coefficient of y^4: "<< get_coefficient(F,4) << std::endl;
  std::cout << "Coefficient of y^5: "<< get_coefficient(F,5) << std::endl; 
  std::cout << std::endl;
  
  PT_2::Leading_coefficient lcoeff;
  std::cout << "Leading coefficient with respect to y:           "
            << lcoeff(F)   // = 11*x^2 + 5*x
            << std::endl;

  PT_2::Get_innermost_coefficient get_icoeff;
  std::cout << "Innermost coefficient of monomial x^1y^4:        "
            << get_icoeff(F,CGAL::Exponent_vector(1,4)) // = 5
            << std::endl;
  
  PT_2::Innermost_leading_coefficient ilcoeff;
  std::cout << "Innermost leading coefficient with respect to y: "
            << ilcoeff(F) // = 11 
            << std::endl; 
}
