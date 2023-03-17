#include <CGAL/_test_fct_ch_I_2.h>

#include <CGAL/Homogeneous.h>

#include <CGAL/convex_hull_constructive_traits_2.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#endif// CGAL_USE_LEDA

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

int main()
{
#ifdef CGAL_USE_LEDA
  CGAL::Convex_hull_constructive_traits_2< CGAL::Homogeneous<leda_integer> > cch_H_integer;
  std::cout << "Homogeneous<integer>:" << std::endl;
  CGAL::ch__batch_test(cch_H_integer);
#endif

#ifdef CGAL_USE_GMP
  CGAL::Convex_hull_constructive_traits_2< CGAL::Homogeneous<CGAL::Gmpz> > cch_H_gmp;
  std::cout << "Homogeneous<gmp>:" << std::endl;
  CGAL::ch__batch_test(cch_H_gmp);
#endif

  return EXIT_SUCCESS;
}
