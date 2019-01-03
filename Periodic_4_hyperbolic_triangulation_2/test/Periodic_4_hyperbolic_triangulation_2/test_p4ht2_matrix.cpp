#include <CGAL/CORE_Expr.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/internal/Exact_complex.h>
#include <CGAL/internal/Hyperbolic_octagon_translation_matrix.h>

#include <iostream>
#include <vector>

using namespace CGAL;
using namespace std;

int main(int, char**)
{
  typedef CORE::Expr                                   NT;
  typedef Exact_complex<NT>                            ECplx;
  typedef Hyperbolic_octagon_translation_matrix<ECplx> Matrix;
  typedef Point_2< Cartesian<NT> >                     Point;

  Matrix m;
  std::cout << "Identity matrix: " << m << std::endl;

  vector<Matrix> gens;
  Matrix::generators(gens);
  for(std::size_t i=0; i<gens.size(); ++i)
    std::cout << "g[" << i << "] = " << gens[i] << std::endl;

  CGAL_assertion(gens[0]*gens[4] == m);
  CGAL_assertion(gens[1]*gens[5] == m);
  CGAL_assertion(gens[2]*gens[6] == m);
  CGAL_assertion(gens[3]*gens[7] == m);

  ECplx o(NT(0), NT(0));
  vector<ECplx> imp;
  for(std::size_t i=0; i<gens.size(); ++i)
  {
    imp.push_back(gens[i](o));
    std::cout << "imp[" << i << "] = " << imp[i] << std::endl;
  }

  CGAL_assertion(imp[0] == ECplx(CGAL::sqrt(NT(2))/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2)))), NT(0)));
  CGAL_assertion(imp[1] == ECplx(NT(1)/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2)))), NT(1)/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2))))));

  std::cout << "test concluded successfully!" << std::endl;

  return EXIT_SUCCESS;
}
