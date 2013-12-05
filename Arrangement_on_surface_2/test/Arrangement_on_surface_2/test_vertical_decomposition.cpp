#include <iostream>
#include <boost/lexical_cast.hpp>

#include <CGAL/basic.h>

#include "test_configuration.h"

#if ((TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS) ||	\
     (TEST_GEOM_TRAITS == BEZIER_GEOM_TRAITS) ||	\
     (TEST_GEOM_TRAITS == RATIONAL_ARC_GEOM_TRAITS)) && !defined(CGAL_USE_CORE)

int main()
{
//  bool   UNTESTED_TRAITS_AS_CORE_IS_NOT_INSTALLED;
  std::cout << std::endl
            << "NOTE: Core is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}
#elif (TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS) && \
  (TEST_NT == LEDA_INT_NT || TEST_NT == LEDA_RAT_NT) && \
  (! CGAL_USE_LEDA)

int main()
{
//  bool   UNTESTED_TRAITS_AS_LEDA_IS_NOT_INSTALLED;
  std::cout << std::endl
	    << "NOTE: LEDA is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}

#elif (TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS) && \
  (TEST_NT == CGAL_GMPZ_NT || TEST_NT == CGAL_GMPQ_NT) && \
  ! (CGAL_USE_GMP && CGAL_USE_MPFI)

int main()
{

//  bool   UNTESTED_TRAITS_AS_GMP_OR_MPFI_IS_NOT_INSTALLED;
  std::cout << std::endl
	    << "NOTE: GMP and/or MPFI are not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}

#elif (TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS) && \
  (TEST_NT == CORE_INT_NT) && \
  !CGAL_USE_CORE

int main()
{
//  bool   UNTESTED_TRAITS_AS_CORE_IS_NOT_INSTALLED;
  std::cout << std::endl
	    << "NOTE: CORE is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}


#else

#include "test_traits.h"
#include "Vertical_decomposition_test.h"

bool test(const char* points_filename, const char* xcurves_filename,
          const char* curves_filename, size_t verbose_level)
{
  Geom_traits geom_traits;
  Vertical_decomposition_test<Geom_traits, Topol_traits> pl_test(geom_traits);
  pl_test.set_verbose_level(verbose_level);
  pl_test.set_filenames(points_filename, xcurves_filename, curves_filename);

  if (!pl_test.init()) return false;
  if (!pl_test.perform()) return false;
  pl_test.clear();

  return true;
}

int main(int argc, char* argv[])
{
#if TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
#endif

  size_t verbose_level = 0;
  int success = 0;
  size_t i = 1;

  // Test 1
  if ((argc > 2) && (std::strncmp(argv[1], "-v", 2) == 0)) {
    verbose_level = boost::lexical_cast<size_t>(argv[2]);
    i += 2;
  }

  if (argc < static_cast<int>(i + 3)) {
    std::cout << "Usage: " << argv[0]
              << " point-file xcurve-file curve-file"
              << std::endl;
    std::cout << "point-file   - the input point file" << std::endl;
    std::cout << "xcurve-file  - the input x-monotone curves file" << std::endl;
    std::cout << "curve-file   - the input curve file" << std::endl;
    return -1;
  }

  for (; static_cast<int>(i) < argc; i += 3) {
    const char* points_filename = argv[i];
    const char* xcurves_filename = argv[i+1];
    const char* curves_filename = argv[i+2];

    if (!test(points_filename, xcurves_filename, curves_filename, verbose_level))
    {
      std::cout << "ERROR : " << argv[0] << " "
                << points_filename << " " << xcurves_filename << " "
                << curves_filename << std::endl;
      success = -1;
    }
  }
  return success;
}

#endif
