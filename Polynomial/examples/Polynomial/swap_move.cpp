#include <CGAL/config.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  typedef CGAL::Polynomial_type_generator<int,3>::Type Poly_3;
  typedef CGAL::Polynomial_traits_d<Poly_3>            PT_3;
  
  //construction using shift 
  Poly_3 x = PT_3::Shift()(Poly_3(1),1,0); // x_0^1
  Poly_3 y = PT_3::Shift()(Poly_3(1),1,1); // x_1^1
  Poly_3 z = PT_3::Shift()(Poly_3(1),1,2); // x_2^1
  
  
  Poly_3 F = 3*x + 5*y + 7*z; 
  std::cout << "The trivariate polynomial F: " << F << std::endl;
  std::cout << std::endl;

  PT_3::Swap swap; 
  PT_3::Move move; 

  std::cout << "x (x_0) and z (x_2) swapped: "<< swap(F,0,2) <<std::endl;  
  std::cout << "x (x_0) and y (x_1) swapped: "<< swap(F,0,1) <<std::endl;
  std::cout << std::endl;
  std::cout << "x (x_0) moved to outermost position (x_2): "<< move(F,0,2) 
            << std::endl;
  std::cout << "Same as swap(swap(F,0,1),1,2): "<< swap(swap(F,0,1),1,2) 
            << std::endl;
  
}
