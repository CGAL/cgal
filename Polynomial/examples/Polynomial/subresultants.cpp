#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Exact_integer.h>

int main(){
  CGAL::set_pretty_mode(std::cout);

  typedef CGAL::Exact_integer Int;

  typedef CGAL::Polynomial_type_generator<Int,1>::Type Poly_1;
  typedef CGAL::Polynomial_traits_d<Poly_1>            PT_1;

  //construction using shift
  Poly_1 x = PT_1::Shift()(Poly_1(1),1); // x^1

  Poly_1 F // = (x+1)^2*(x-1)*(2x-1)=2x^4+x^3-3x^2-x+1
    =   2 * CGAL::ipower(x,4) + 1 * CGAL::ipower(x,3)
      - 3 * CGAL::ipower(x,2) - 1 * CGAL::ipower(x,1)
      + 1 * CGAL::ipower(x,0);
  std::cout << "F=" << F << std::endl;

  Poly_1 G // = (x+1)*(x+3)=x^2+4*x+3
    =   1 * CGAL::ipower(x,2) + 4 * CGAL::ipower(x,1) + 3 * CGAL::ipower(x,0);
  std::cout << "G=" << G << std::endl;

  // Resultant computation:
  PT_1::Resultant resultant;

  std::cout << "The resultant of F and G is: " << resultant(F,G) << std::endl;
  // It is zero, because F and G have a common factor

  // Real root counting:
  PT_1::Principal_sturm_habicht_sequence stha;
  std::vector<Int> psc;

  stha(F,std::back_inserter(psc));

  int roots = CGAL::number_of_real_roots(psc.begin(),psc.end());

  std::cout << "The number of real roots of F is: " << roots << std::endl; // 3

  roots =  CGAL::number_of_real_roots(G);

  std::cout << "The number of real roots of G is: " << roots << std::endl; // 2

  return 0;

}
