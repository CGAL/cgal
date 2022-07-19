#include <CGAL/_test_fct_ch_I_2.h>

#include <CGAL/Homogeneous.h>

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
  CGAL::Homogeneous<leda_integer> ch_H_integer;
  std::cout << "Homogeneous<integer>:" << std::endl;
  CGAL::ch__batch_test( ch_H_integer );
#endif

#ifdef CGAL_USE_GMP
  CGAL::Homogeneous<CGAL::Gmpz> ch_H_gmp;
  std::cout << "Homogeneous<gmp>:" << std::endl;
  CGAL::ch__batch_test( ch_H_gmp );
#endif

  CGAL::Homogeneous<double> ch_H_double;
  std::cout << "Homogeneous<double>:" << std::endl;
  CGAL::ch__batch_test( ch_H_double );

  return EXIT_SUCCESS;
}
