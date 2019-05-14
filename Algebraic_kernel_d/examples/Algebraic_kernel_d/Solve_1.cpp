// $URL$
// $Id$

#include <CGAL/basic.h>
#ifdef CGAL_USE_MPFI 
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Gmpz.h>
#include <vector>

typedef CGAL::Algebraic_kernel_d_1<CGAL::Gmpz>          AK;
typedef AK::Polynomial_1                                Polynomial_1;
typedef AK::Algebraic_real_1                            Algebraic_real_1;
typedef AK::Bound                                       Bound;
typedef AK::Multiplicity_type                           Multiplicity_type;

int main(){
  AK ak; // an object of 
  AK::Solve_1 solve_1 = ak.solve_1_object();
  Polynomial_1 x = CGAL::shift(AK::Polynomial_1(1),1); // the monomial x


  // variant using a bool indicating a square free polynomial
  // multiplicities are not computed
  std::vector<Algebraic_real_1> roots;
  solve_1(x*x-2,true, std::back_inserter(roots));
  std::cout << "Number of roots is           : " << roots.size() << "\n";
  std::cout << "First root should be -sqrt(2): " << CGAL::to_double(roots[0]) << "\n";
  std::cout << "Second root should be sqrt(2): " << CGAL::to_double(roots[1]) << "\n\n";
  roots.clear();

  // variant for roots in a given range of a square free polynomial
  solve_1((x*x-2)*(x*x-3),true, Bound(0),Bound(10),std::back_inserter(roots));
  std::cout << "Number of roots is           : " << roots.size() << "\n";
  std::cout << "First root should be  sqrt(2): " << CGAL::to_double(roots[0]) << "\n";
  std::cout << "Second root should be sqrt(3): " << CGAL::to_double(roots[1]) << "\n\n";
  roots.clear();

  // variant computing all roots with multiplicities
  std::vector<std::pair<Algebraic_real_1,Multiplicity_type> > mroots;
  solve_1((x*x-2), std::back_inserter(mroots));
  std::cout << "Number of roots is           : " << mroots.size() << "\n";
  std::cout << "First root should be -sqrt(2): " << CGAL::to_double(mroots[0].first) << ""
            << " with multiplicity "             << mroots[0].second << "\n";
  std::cout << "Second root should be sqrt(2): " << CGAL::to_double(mroots[1].first) << ""
            << " with multiplicity "             << mroots[1].second << "\n\n";
  mroots.clear();

  // variant computing roots with multiplicities for a range
  solve_1((x*x-2)*(x*x-3),Bound(0),Bound(10),std::back_inserter(mroots));
  std::cout << "Number of roots is           : " << mroots.size() << "\n";
  std::cout << "First root should be  sqrt(2): " << CGAL::to_double(mroots[0].first) << ""
            << " with multiplicity "             << mroots[0].second << "\n";
  std::cout << "Second root should be sqrt(3): " << CGAL::to_double(mroots[1].first) << ""
            << " with multiplicity "             << mroots[1].second << "\n\n";
  return 0;
}
#else
int main(){
  std::cout << "This example requires CGAL to be configured with library MPFI." << std::endl;
return 0;
}
#endif
