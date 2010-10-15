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
  
  
  Poly_3 F = x*y*y*z*z*z; 
  std::cout << "The trivariate polynomial F: " << F << std::endl;
  std::cout << std::endl;

  PT_3::Swap swap; 
  PT_3::Move move; 
  PT_3::Permute permute; 

  std::cout << "x and z swapped: "<< swap(F,0,2) // = x^3*y^2*z
            << std::endl;  
  std::cout << "x and y swapped: "<< swap(F,0,1) // = x^2*y*z^3
            << std::endl << std::endl; 
  
  std::cout << "x moved to outermost position           : "
            << move(F,0,2)                       // = x^2*y^3*z
            << std::endl;
  std::cout << "Same as swap(swap(F,0,1),1,2)           : "
            << swap(swap(F,0,1),1,2)             // = x^2*y^3*z
            << std::endl;

  std::cout << "Same as the permutation (0,1,2)->(2,0,1): ";
  std::vector<int> perm; 
  perm.push_back(2);perm.push_back(0);perm.push_back(1);
  std::cout << permute(F,perm.begin(),perm.end())// = x^2*y^3*z
            << std::endl;
  
}
