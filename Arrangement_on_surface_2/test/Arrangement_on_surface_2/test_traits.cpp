#include <CGAL/basic.h>
#include "test_configuration.h"
#include <iostream>

#if ((TEST_TRAITS == CORE_CONIC_TRAITS) ||	\
     (TEST_TRAITS == BEZIER_TRAITS) ||					\
     (TEST_TRAITS == RATIONAL_ARC_TRAITS)) && !defined(CGAL_USE_CORE)

int main ()
{
  bool   UNTESTED_TRAITS_AS_CORE_IS_NOT_INSTALLED;
  std::cout << std::endl
            << "WARNING: Core is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}
#elif (TEST_TRAITS == ALGEBRAIC_TRAITS) && \
  (TEST_NT == LEDA_INT_NT || TEST_NT == LEDA_RAT_NT) && \
  (! CGAL_USE_LEDA)

int main ()
{
  bool   UNTESTED_TRAITS_AS_LEDA_IS_NOT_INSTALLED;
  std::cout << std::endl
	    << "LEDA is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}

#elif (TEST_TRAITS == ALGEBRAIC_TRAITS) && \
  (TEST_NT == CGAL_GMPZ_NT || TEST_NT == CGAL_GMPQ_NT) && \
  ! (CGAL_USE_GMP && CGAL_USE_MPFI)

int main ()
{

  bool   UNTESTED_TRAITS_AS_GMP_OR_MPFI_IS_NOT_INSTALLED;
  std::cout << std::endl
	    << "GMP and/or MPFI are not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}

#elif (TEST_TRAITS == ALGEBRAIC_TRAITS) && \
  (TEST_NT == CORE_INT_NT) && \
  !CGAL_USE_CORE

int main ()
{
  bool   UNTESTED_TRAITS_AS_CORE_IS_NOT_INSTALLED;
  std::cout << std::endl
	    << "CORE is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}


#else

#include <CGAL/assertions.h>
#include <CGAL/Arrangement_2.h>

#include "test_traits.h"
#include "Traits_test.h"

int main (int argc, char * argv[])
{
#if TEST_TRAITS == ALGEBRAIC_TRAITS
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
#endif

  Traits_test<Traits> test(argc, argv);
  bool rc  = test.start();
  return (rc) ? 0 : -1;
}

#endif
