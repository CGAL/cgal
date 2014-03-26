// #define CGAL_OVERLAY_TRAITS_VERBOSE 1

/* Test the overlay operation.
 * 1. Read a pair of sets---a set of x-monotone-curves followed by a set of
 *    (isolated) vertices.
 * 2. Read another pair of such sets.
 * 3. Read the expected extended arrangement and store in m_arr.
 *    Each type of cell, i.e., vertex, halfedge, and face, is extended with
 *    an unsigned int.
 * 4. Construct two extended arrangements, namely m_arr1 and m_arr2, induced
 *    by the elements in the two pairs above, respectively.
 * 5. Initialize the data field of each halfedge with the number of curves
 *    that induced that halfedge if the halfedge is directed left-to-right
 *    and twice the the number of curves that induced that halfedge if the
 *    halfedge is directed right-to-left. We initialize the data field of
 *    each face with the total sum of the data of the halfedges on the
 *    boundary of the face. We initialize the data field of each isolated
 *    vertex with the data of the containing face, and the data field of each
 *    non-isolated vertex with the total sum of the data of the incident
 *    halfedges (the vertex is their target.
 * 6. Compute the overlay of arr1 and arr2 and store the result in a third
 *    arrangement called arr. The data field of each cell of arr is assigned
 *    with the sum of the data of the 2 inducing cells from arr1 and arr2,
 *    respectively.
 * 7. Finally, test whether arr and m_arr are equivalent (isomorphic).
 */

#include <iostream>
#include <cstring> // for std::strncmp
#include <string.h>
#include <stdlib.h>

#include "test_configuration.h"

#include <CGAL/basic.h>

#if ((TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS) ||	\
     (TEST_GEOM_TRAITS == BEZIER_GEOM_TRAITS) ||	\
     (TEST_GEOM_TRAITS == RATIONAL_ARC_GEOM_TRAITS)) && !defined(CGAL_USE_CORE)

int main()
{
  // bool UNTESTED_TRAITS_AS_CORE_IS_NOT_INSTALLED;
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
  // bool UNTESTED_TRAITS_AS_LEDA_IS_NOT_INSTALLED;
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

  // bool UNTESTED_TRAITS_AS_GMP_OR_MPFI_IS_NOT_INSTALLED;
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
  //  bool UNTESTED_TRAITS_AS_CORE_IS_NOT_INSTALLED;
  std::cout << std::endl
	    << "NOTE: CORE is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}

#else

#include "test_geom_traits.h"

#include <CGAL/Arr_curve_data_traits_2.h>

// Define Geom_traits to be the curve-data-traits of the base geom traits.
typedef CGAL::Arr_curve_data_traits_2<Base_geom_traits,
                                      unsigned int,
                                      std::plus<unsigned int> >  
                                                        Geom_traits;
typedef Geom_traits::Point_2                            Point_2;
typedef Geom_traits::Curve_2                            Curve_2;
typedef Geom_traits::X_monotone_curve_2                 X_monotone_curve_2;

#include <CGAL/Arr_extended_dcel.h>

typedef CGAL::Arr_extended_dcel<Geom_traits,
                                unsigned int, unsigned int, unsigned int>
                                                          Dcel;

#include "test_extended_topol_traits.h"
#include "Overlay_test.h"

bool test(const char* filename, int verbose_level)
{
  Geom_traits geom_traits;
  Overlay_test<Geom_traits, Topol_traits> overlay_test(geom_traits);
  overlay_test.set_verbose_level(verbose_level);
  overlay_test.set_filename(filename);
  if (!overlay_test.init()) return false;
  if (!overlay_test.perform()) return false;
  overlay_test.clear();
  return true;
}

int main(int argc, char* argv[])
{
#if TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
#endif

  if (argc < 2) {
    std::cerr << "Missing input file!" << std::endl;
    std::cerr << "Usage: " << argv[0] << " input file" << std::endl;
    return -1;
  }

  // TBD: Replace with better parsing!
  int i = 1;
  int verbose_level = 0;
  if (argc > 2) {
    if ((argc > 3) && (std::strncmp(argv[1], "-v", 2) == 0)) {
      verbose_level = atoi(argv[2]);
      i += 2;
    }
  }
  
  int success = 0;
  for (; i < argc; ++i) {
    const char* filename = argv[i];

    if (!test(filename, verbose_level)) {
      std::cout << "ERROR : " << argv[0] << " " << filename << std::endl;
      success = -1;
    }
  }
    
  return success;
}

#endif
