#include <CGAL/_test_fct_ch_I_2.h>

#include <CGAL/Cartesian.h>

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
  CGAL::Cartesian<leda_rational> ch_C_rational;
  std::cout << "Cartesian<rational>:" << std::endl;
  CGAL::ch__batch_test(ch_C_rational);
#endif

#ifdef CGAL_USE_GMP
  CGAL::Cartesian<CGAL::Quotient<CGAL::Gmpz> > ch_C_Qgmp;
  std::cout << "Cartesian<Quotient<Gmpz> > >:" << std::endl;
  CGAL::ch__batch_test(ch_C_Qgmp);
#endif

  return EXIT_SUCCESS;
}
