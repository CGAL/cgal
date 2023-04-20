#include <CGAL/_test_fct_ch_I_2.h>

#include <CGAL/Simple_cartesian.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

int main()
{
#ifdef CGAL_USE_LEDA
  CGAL::Simple_cartesian<leda_rational> ch_S_rational;
  std::cout << "SimpleCartesian<rational>:" << std::endl;
  CGAL::ch__batch_test( ch_S_rational );
#endif

#ifdef CGAL_USE_GMP
  CGAL::Simple_cartesian<CGAL::Quotient<CGAL::Gmpz> > ch_S_Qgmp;
  std::cout << "SimpleCartesian<Quotient<Gmpz> > >:" << std::endl;
  CGAL::ch__batch_test( ch_S_Qgmp );
#endif

  return EXIT_SUCCESS;
}
