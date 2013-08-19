#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  typedef CGAL::Polynomial_type_generator<int,2>::Type Poly_2;
  typedef CGAL::Polynomial_traits_d<Poly_2>            PT_2;
  
  //construction using shift 
  Poly_2 x = PT_2::Shift()(Poly_2(1),1,0); // x_0^1
  Poly_2 y = PT_2::Shift()(Poly_2(1),1,1); // x_1^1
  
  Poly_2 F // = (11*x^2 + 5*x)*y^4 + (7*x^2)*y^3
    = 11 * CGAL::ipower(y,4) * CGAL::ipower(x,2) 
    + 5 * CGAL::ipower(y,4)  * CGAL::ipower(x,1) 
    + 7 * CGAL::ipower(y,3)  * CGAL::ipower(x,2);  
  std::cout << "The bivariate polynomial F: " << F <<"\n"<< std::endl;
  
  PT_2::Degree degree; 
  PT_2::Total_degree total_degree; 
  PT_2::Degree_vector degree_vector;
  
  std::cout << "The degree of F with respect to y: "<< degree(F)       // = 4 
            << std::endl;
  std::cout << "The degree of F with respect to x: "<< degree(F,0)     // = 2 
            << std::endl;
  std::cout << "The total degree of F            : "<< total_degree(F) // = 6 
            << std::endl;
  std::cout << "The degree vector of F           : "<< degree_vector(F)// = (2,4)
            << std::endl;
  
  
  
}
