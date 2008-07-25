#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Gmpz.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  typedef CGAL::Polynomial_type_generator<CGAL::Gmpz,3>::Type Poly_int_3;
  typedef CGAL::Polynomial_traits_d<Poly_int_3>            PT_3;
  typedef PT_3::Coefficient                                Poly_int_1;
  typedef PT_3::Innermost_coefficient                      Integer; 
  
  //construction using shift 
  Poly_int_3 x = PT_3::Shift()(Poly_int_3(1),1,0); // shift 1 by x_0^1
  Poly_int_3 y = PT_3::Shift()(Poly_int_3(1),1,1); // shift 1 by x_1^1
  Poly_int_3 z = PT_3::Shift()(Poly_int_3(1),1,2); // shift 1 by x_2^1
  
  
  Poly_int_3 F = 3*x + 5*y + 7*z; 
  std::cout << "The trivariate polynomial F: " << F << std::endl;
  std::cout << std::endl;

  PT_3::Swap swap; 
  PT_3::Move move; 

  std::cout << "x and z swapped: "<< swap(F,0,2) <<std::endl;  
  std::cout << "x and y swapped: "<< swap(F,0,1) <<std::endl;
  std::cout << std::endl;
  std::cout << "x moved to outermost position: "<< move(F,0,2) 
            <<std::endl;
  std::cout << "Same as swap(swap(F,0,1),1,2): "<< swap(swap(F,0,1),1,2) 
            <<std::endl;
  
}
