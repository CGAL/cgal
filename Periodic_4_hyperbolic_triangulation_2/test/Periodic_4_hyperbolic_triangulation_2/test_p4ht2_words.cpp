#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

// #include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/internal/Exact_complex.h>
#include <CGAL/Hyperbolic_octagon_translation.h>

#include <iostream>
#include <vector>

int main(int, char**)
{
  typedef CORE::Expr                                          NT;
  typedef CGAL::Exact_complex<NT>                             ECplx;
  typedef CGAL::Hyperbolic_octagon_translation<NT>            Word;

  Word w;
  // std::cout << "empty word: " << w << ", matrix: " << w.matrix() << std::endl;

  Word a(0), ab(0, 5), abc(0, 5, 2), abcd(0, 5, 2, 7), dcb(7, 2, 5), dc(7, 2), d(7);
  // std::cout << "a    = " << a << ",    matrix: " << a.matrix()     << std::endl;
  // std::cout << "ab   = " << ab << ",   matrix: " << ab.matrix()    << std::endl;
  // std::cout << "abc  = " << abc << ",  matrix: " << abc.matrix()   << std::endl;
  // std::cout << "abcd = " << abcd << ", matrix: " << abcd.matrix()  << std::endl;
  // std::cout << "dcb  = " << dcb << ",  matrix: " << dcb.matrix()   << std::endl;
  // std::cout << "dc   = " << dc << ",   matrix: " << dc.matrix()    << std::endl;
  // std::cout << "d    = " << d << ",    matrix: " << d.matrix()     << std::endl;

  return EXIT_SUCCESS;
}
