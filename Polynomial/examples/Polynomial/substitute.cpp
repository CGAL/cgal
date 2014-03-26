#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  typedef CGAL::Polynomial_type_generator<int,2>::Type Poly_2;
  typedef CGAL::Polynomial_traits_d<Poly_2>            PT_2;
  
  //construction using shift 
  Poly_2 x = PT_2::Shift()(Poly_2(1),1,0); // x^1
  Poly_2 y = PT_2::Shift()(Poly_2(1),1,1); // y^1
  
  
  Poly_2 F = 2*x*y + 3*CGAL::ipower(y,3);
  std::cout << "The bivariate polynomial F: " << F // = 3*y^3 + (2*x)*y
            << std::endl << std::endl;

  PT_2::Evaluate evaluate;
  PT_2::Evaluate_homogeneous hevaluate;
  
  // Evaluation considers a polynomials as univariate:
  std::cout << "F(5): " << evaluate(F,5)      // = 10*x + 375
            << std::endl; 
  // Evaluate_homogeneous considers F as a homogeneous polynomial in 
  // the outermost variable only, that is, F is interpreted as 
  // F(u,v) = 2*x*u*v^2 + 3 * u^3 
  std::cout << "F(5,7): " << hevaluate(F,5,7) // = 490*x + 375
            << std::endl << std::endl;

  
  PT_2::Substitute substitute;
  PT_2::Substitute_homogeneous hsubstitute;
  
  // Substitute considers a polynomials as multivariate, that is, the 
  // new values for the variables are given by an iterator range
  // Note that the value type only has to be interoperable with the innermost 
  // coefficient
  std::list<Poly_2> replacements;
  replacements.push_back(x-1); // replace x by x-1
  replacements.push_back(y);   // replace y by y, i.e., do nothing

  std::cout << "The bivariate polynomial F: " << F // = 3*y^3 + (2*x)*y
            << std::endl;
  std::cout << "F(x-1,y):   " // = 3*y^3 + (2*x + (-2))*y
            << substitute(F,replacements.begin(),replacements.end())
            << std::endl;
  // Substitute_homogeneous considers F as a homogeneous polynomial in 
  // all variable, that is, F is interpreted as 
  // F(x,y,w) = 2*x*y*w + 3 * y^3 
  replacements.push_back(y);  // replace z by y 

  std::cout << "F(x-1,y,y): " // = 3*y^3 + (2*x + (-2))*y^2
            << hsubstitute(F,replacements.begin(),replacements.end())
            << std::endl;
}
