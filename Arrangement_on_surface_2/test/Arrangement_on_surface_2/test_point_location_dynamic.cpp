#include <iostream>

#include "test_configuration.h"

#include <CGAL/config.h>

#if ((TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS) ||        \
     (TEST_GEOM_TRAITS == BEZIER_GEOM_TRAITS) ||        \
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
#include "Point_location_dynamic_test.h"

bool test1(const char* points_filename, const char* xcurves_filename,
           const char* curves_filename, const char* commands_filename,
           const char* queries_filename)
{
  Geom_traits geom_traits;
  Point_location_dynamic_test<Geom_traits, Topol_traits> pl_test(geom_traits);
  pl_test.set_filenames(points_filename, xcurves_filename, curves_filename,
                        commands_filename, queries_filename);

  if (!pl_test.allocate_arrangement()) return false;
  if (!pl_test.construct_pl_strategies()) return false;
  if (!pl_test.init()) return false;
  if (!pl_test.construct_arrangement()) return false;
  if (!pl_test.perform()) return false;

  pl_test.clear();
  pl_test.deallocate_arrangement();
  pl_test.deallocate_pl_strategies();

  return true;
}

int main(int argc, char* argv[])
{
#if TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS
  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);
#endif

  if (argc < 5) {
    std::cout << "Usage: " << argv[0]
              << " point-file xcurve-file curve-file command-file query-file"
              << std::endl;
    std::cout << "point-file   - the input point file" << std::endl;
    std::cout << "xcurve-file  - the input x-monotone curves file" << std::endl;
    std::cout << "command-file - the command file" << std::endl;
    std::cout << "curve-file   - the input curve file" << std::endl;
    std::cerr << "query-file  - the input query point file" << std::endl;
    return -1;
  }

  int success = 0;

  // Test 1
  for (int i = 1; i < argc; i += 5) {
    const char* points_filename = argv[i];
    const char* xcurves_filename = argv[i+1];
    const char* curves_filename = argv[i+2];
    const char* commands_filename = argv[i+3];
    const char* queries_filename = argv[i+4];

    if (!test1(points_filename, xcurves_filename,
               curves_filename, commands_filename, queries_filename))
    {
      std::cout << "ERROR : " << argv[0] << " " << points_filename << " "
                << xcurves_filename << " " << curves_filename << " "
                << commands_filename << " " << queries_filename << std::endl;
      success = -1;
    }
  }
  return success;
}

#endif
