#include <CGAL/basic.h>
#include "test_configuration.h"
#include <iostream>

#if ((TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS) || \
     (TEST_GEOM_TRAITS == BEZIER_GEOM_TRAITS) || \
     (TEST_GEOM_TRAITS == RATIONAL_ARC_GEOM_TRAITS)) && !defined(CGAL_USE_CORE)

int main()
{
//  bool   UNTESTED_GEOM_TRAITS_AS_CORE_IS_NOT_INSTALLED;
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
//  bool   UNTESTED_GEOM_TRAITS_AS_LEDA_IS_NOT_INSTALLED;
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

//  bool   UNTESTED_GEOM_TRAITS_AS_GMP_OR_MPFI_IS_NOT_INSTALLED;
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
//  bool   UNTESTED_GEOM_TRAITS_AS_CORE_IS_NOT_INSTALLED;
  std::cout << std::endl
	    << "NOTE: CORE is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}


#else

#include "test_geom_traits.h"
// Define Base_geom_traits to be the geometry traits, namely, Geom_traits.
typedef Base_geom_traits                Geom_traits;
typedef Geom_traits::Point_2            Point_2;
typedef Geom_traits::Curve_2            Curve_2;
typedef Geom_traits::X_monotone_curve_2 X_monotone_curve_2;
#include "Traits_test.h"

int main(int argc, char* argv[])
{

#if TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
#endif

  Geom_traits traits;

  Traits_test<Geom_traits> test(traits);

  /* the parse() function is present in "Traits_base_tests.h" which further
   * call parse() in "IO_test.h"
   * The parse() in "Traits_base_test.h" simply assign the command file which
   * is argv[4] i.e. compare_y_at_x file to the 'm_filename_commands' variable
   * and also assigns the polycurve_conic_traits string.
   * The parse() in "IO_test.h" assigns the argv[1, 2 and 3] i.e files paths
   * of points, xcurves and curves to appropriate string variables.
   *
   * Inheritance tree
   *   Traits_test <-- Traits_base_test <-- IO_test <-- IO_base_test
   */
  if (!test.parse(argc, argv)) return -1;

  /* Inheritance tree
   *   Traits_test <-- Traits_base_test <-- IO_test <-- IO_base_test
   *
   * The init() function is present in "IO_test.h" the function simply initiate
   * the points, curves and x-monotone curves vectors through the input data
   * files provided by calling read_points() read_curves() and read_xcurves()
   * respectively. Each of these functions calls read_point(), read_curve()
   * and read_xcurve() for construction from "IO_base_test.h"
   *
   * read_point(), read_curve() and read_xcurve() from "IO_base_test.h"
   * construct the appriopriate point, curve and xcurve using appropriate
   * GEOM_TRAITS i.e. in this case POLYCURVE_CONIC_GEOM_TRAITS and using
   * overridden functions. Note: these functions only make 1 curve. So if we
   * want a polycurve, it should be taken care of in these function.
  */
  if (!test.init()) return -1;  //polycurves are made here.

  /* The perform() function is present in Traits_base_test.h
   */
  if (!test.perform()) return -1;

  return 0;
}

#endif
