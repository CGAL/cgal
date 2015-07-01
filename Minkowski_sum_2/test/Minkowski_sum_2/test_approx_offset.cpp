#include <CGAL/Exact_rational.h>

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational         Rational;

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/Polygon_convex_decomposition_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <cstring> // for std::strchr

#include "read_polygon.h"

typedef CGAL::Cartesian<Rational>                 Kernel;
typedef Kernel::Point_2                           Point_2;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;

typedef CGAL::Gps_circle_segment_traits_2<Kernel> Gps_traits_2;
typedef Gps_traits_2::Polygon_2                   Offset_polygon_2;
typedef Gps_traits_2::Polygon_with_holes_2        Offset_polygon_with_holes_2;

/*! Check if two polygons with holes are the same. */
bool are_equal (const Offset_polygon_with_holes_2& ph1,
                const Offset_polygon_with_holes_2& ph2)
{
  std::list<Offset_polygon_with_holes_2>   sym_diff;
  CGAL::symmetric_difference (ph1, ph2, std::back_inserter(sym_diff));
  return (sym_diff.empty());
}

/*! The main program. */
int main (int argc, char* argv[])
{
  const double eps = 0.0001;

  // Read the input file. Because of the structure of the *.cmd file 
  // (which is concatenated to the command line) we need to get all the
  // inputs in one command line. This is the reason we read triplets/pairs of
  // arguments. Each triplet/double is one input for the program.
  if (argc < 3) {
    std::cerr << "Usage (input is triplets of): <polygon> <radius> "
      "[decomposition flags]" << std::endl;
    return -1;
  }

  int i = 1;
  while (i < argc) {
    // Read the polygon from the input file.
    Polygon_2   pgn;
    const char* filename = argv[i];
    read_polygon (filename, pgn);
    
    // Read the offset radius.
    Rational r;
    std::istringstream iss (argv[i+1], std::istringstream::in);
    iss >> r;

    std::cout << "Testing " << filename << " and " << r << std::endl;

    // Read the decomposition flags.
    bool use_ssab = true;
    bool use_opt = true;
    bool use_hm = true;
    bool use_greene = true;
    
    if (((i+2) < argc) && (argv[i+2][0] == '-')) {
      use_ssab = (std::strchr (argv[i+2], 's') != NULL);
      use_opt = (std::strchr (argv[i+2], 'o') != NULL);
      use_hm = (std::strchr (argv[i+2], 'h') != NULL);
      use_greene = (std::strchr (argv[i+2], 'g') != NULL);
    }
    
    // Compute the Minkowski sum using the convolution method.
    Offset_polygon_with_holes_2 offset_conv;
    
    std::cout << "Using the convolution method ... ";
    offset_conv = approximated_offset_2 (pgn, r, eps);
    std::cout << "Done." << std::endl;
    
    // Define auxiliary polygon-decomposition objects.
    CGAL::Small_side_angle_bisector_decomposition_2<Kernel> ssab_decomp;
    CGAL::Optimal_convex_decomposition_2<Kernel>            opt_decomp;
    CGAL::Hertel_Mehlhorn_convex_decomposition_2<Kernel>    hm_approx_decomp;
    CGAL::Greene_convex_decomposition_2<Kernel>             greene_decomp;
    Offset_polygon_with_holes_2                             offset_decomp;
    
    if (use_ssab) {
      std::cout << "Using the small-side angle-bisector decomposition ... ";
      offset_decomp = approximated_offset_2 (pgn, r, eps, ssab_decomp);
      if (are_equal (offset_conv, offset_decomp))
        std::cout << "OK." << std::endl;
      else
        std::cout << "ERROR (different result)." << std::endl;
    }
    
    if (use_opt) {
      std::cout << "Using the optimal convex decomposition ... ";
      offset_decomp = approximated_offset_2 (pgn, r, eps, opt_decomp);
      if (are_equal (offset_conv, offset_decomp))
        std::cout << "OK." << std::endl;
      else
        std::cout << "ERROR (different result)." << std::endl;
    }
    
    if (use_hm) {
      std::cout << "Using the Hertel--Mehlhorn decomposition ... ";
      offset_decomp = approximated_offset_2 (pgn, r, eps, hm_approx_decomp);
      if (are_equal (offset_conv, offset_decomp))
        std::cout << "OK." << std::endl;
      else
        std::cout << "ERROR (different result)." << std::endl;
    }
    
    if (use_greene) {
      std::cout << "Using the Greene decomposition ... ";
      offset_decomp = approximated_offset_2 (pgn, r, eps, greene_decomp);
      if (are_equal (offset_conv, offset_decomp))
        std::cout << "OK." << std::endl;
      else
        std::cout << "ERROR (different result)." << std::endl;
    }
    
    i += (((i+2) < argc) && (argv[i+2][0] == '-')) ? 3 : 2;
  }
  
  return 0;
}
