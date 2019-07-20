#include <CGAL/Exact_algebraic.h>
#include <CGAL/internal/Exact_complex.h>
#include <CGAL/internal/Hyperbolic_octagon_translation_matrix.h>
#include <iostream>
#include <vector>

int main(int, char**)
{

  typedef CGAL::Exact_algebraic                               NT;
  typedef CGAL::Exact_complex<NT>                             ECplx;
  typedef CGAL::Hyperbolic_octagon_translation_matrix<ECplx>  Matrix;

  Matrix m;
  std::cout << "Identity matrix: " << m << std::endl;

  std::vector<Matrix> gens;
  Matrix::generators(gens);
  for(std::size_t i=0; i<gens.size(); ++i)
    std::cout << "g[" << i << "] = " << gens[i] << std::endl;

  assert(gens[0]*gens[4] == m);
  assert(gens[1]*gens[5] == m);
  assert(gens[2]*gens[6] == m);
  assert(gens[3]*gens[7] == m);

  ECplx o(NT(0), NT(0));
  std::vector<ECplx> imp;
  for(std::size_t i=0; i<gens.size(); ++i)
  {
    imp.push_back(gens[i](o));
    std::cout << "imp[" << i << "] = " << imp[i] << std::endl;
  }

  assert(imp[0] == ECplx(CGAL::sqrt(NT(2))/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2)))), NT(0)));
  assert(imp[1] == ECplx(NT(1)/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2)))), NT(1)/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2))))));

  std::cout << "test concluded successfully!" << std::endl;

  return EXIT_SUCCESS;
}
