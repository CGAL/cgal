
#include <CGAL/config.h>
#include <iostream>

// That one should not be needed in the long term:
#include <CGAL/double.h>
#include <CGAL/float.h>
#include <CGAL/int.h>
#include <CGAL/Interval_nt.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif


#include <CGAL/float.h>
#include <CGAL/double.h>
#include <CGAL/int.h>

#if 0
#ifdef CGAL_USE_CLN
#include <CGAL/CLN/cl_integer.h>
// ...
#endif // CLN
#endif

// Should be tested with different current rounding mode ?

template <typename NT>
void
test(const NT &)
{
  NT zero (0);
  NT one (1);
  NT one_third = one / NT(3);

  std::pair<double,double> zero_i      = CGAL_NTS to_interval(zero);
  std::pair<double,double> one_i       = CGAL_NTS to_interval(one);
  std::pair<double,double> one_third_i = CGAL_NTS to_interval(one_third);

  if (zero_i.first > 0 || zero_i.second < 0)
    std::cout << "  BUG zero ! : " << zero_i.first << " " << zero_i.second << std::endl;
  if (one_i.first > 1 || one_i.second < 1)
    std::cout << "  BUG one ! : " << one_i.first << " " << one_i.second << std::endl;
  std::cout << one_third_i.first << " " << one_third_i.second
            << " (not correct for integer types)" << std::endl;
}

int main()
{
  std::cout << "Test program for the to_interval() function." << std::endl;
  std::cout.precision(20);

  std::cout << "Testing double :" << std::endl;
  ::test(double());
  std::cout << "Testing float :" << std::endl;
  ::test(float());
  std::cout << "Testing int :" << std::endl;
  ::test(int());

#ifdef CGAL_USE_LEDA
  std::cout << "Testing leda_real :" << std::endl;
  ::test(leda_real());
  std::cout << "Testing leda_bigfloat :" << std::endl;
  ::test(leda_bigfloat());
  std::cout << "Testing leda_integer :" << std::endl;
  ::test(leda_integer());
  std::cout << "Testing leda_rational :" << std::endl;
  ::test(leda_rational());
#endif

#ifdef CGAL_USE_GMP
  std::cout << "Testing Gmpz :" << std::endl;
  ::test(CGAL::Gmpz());
#endif

#if 0
#ifdef CGAL_USE_CLN
  std::cout << "Testing cl_integer :" << std::endl;
  ::test(cl_integer());
#endif
#endif

  return 0;
}
