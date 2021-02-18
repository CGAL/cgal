#include <CGAL/config.h>

#ifndef CGAL_USE_CORE

#include <iostream>
#include <sstream>

int main ()
{
  std::cout << "Sorry, this test needs CORE ..." << std::endl;
  return (0);
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/offset_polygon_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/Polygon_convex_decomposition_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include "read_polygon.h"
#include <list>


typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;

// instead of
//typedef CGAL::Cartesian<Rational>                       Rat_kernel;
//typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
//typedef CGAL::Arr_conic_traits_2<Rat_kernel,
//                                 Alg_kernel,Nt_traits>  Conic_traits_2;
// workaround for VC++
struct Rat_kernel : public CGAL::Cartesian<Rational> {};
struct Alg_kernel : public CGAL::Cartesian<Algebraic> {};
struct Conic_traits_2 : public CGAL::Arr_conic_traits_2<Rat_kernel,
      Alg_kernel,Nt_traits> {};

typedef Rat_kernel::Point_2                             Point_2;
typedef CGAL::Polygon_2<Rat_kernel>                     Polygon_2;

typedef CGAL::Gps_traits_2<Conic_traits_2>              Gps_traits_2;
typedef Gps_traits_2::Polygon_2                         Offset_polygon_2;
typedef Gps_traits_2::Polygon_with_holes_2   Offset_polygon_with_holes_2;

/*! Check if two polygons with holes are the same. */
bool are_equal (const Offset_polygon_with_holes_2& ph1,
                const Offset_polygon_with_holes_2& ph2)
{
  std::list<Offset_polygon_with_holes_2>   sym_diff;

  CGAL::symmetric_difference (ph1, ph2,
                              std::back_inserter(sym_diff));

  return (sym_diff.empty());
}

/*! The main program. */
int main (int argc, char **argv)
{
  // Read the input file. Because of the structure of the *.cmd file
  // (which is concatenated to the command line) we need to get all the
  // inputs in one command line. This is the reason we read triplets/pairs of
  // arguments. Each triplet/double is one input for the program.
  if (argc < 3)
  {
    std::cerr << "Usage (input is triplets of): <polygon> <radius> "
      "[decomposition flags]" << std::endl;
    return (1);
  }

  // Read the polygon from the input file.
  Polygon_2   pgn;

  int i = 1;
  while (i < argc)
  {
    read_polygon (argv[i], pgn);

    // Read the offset radius.
    int         numer=0, denom = 0;
    std::istringstream iss(argv[i+1]);
    char c;
    iss >> numer >> c >> denom;

    if (iss.bad())
    {
      std::cerr << "Invalid radius: " << argv[i+1] << std::endl;
      return (1);
    }

    std::cout << "Testing " << argv[i] << " with radius " << argv[i+1] <<
      std::endl;

    Rational    r = Rational (numer, denom);

    // Read the decomposition flags.
    bool         use_ssab = true;
    bool         use_opt = true;
    bool         use_hm = true;
    bool         use_greene = true;

    if (i+2 < argc && argv[i+2][0] == '-')
    {
      use_ssab = (strchr (argv[i+2], 's') != nullptr);
      use_opt = (strchr (argv[i+2], 'o') != nullptr);
      use_hm = (strchr (argv[i+2], 'h') != nullptr);
      use_greene = (strchr (argv[i+2], 'g') != nullptr);
    }

    // Compute the Minkowski sum using the convolution method.
    Conic_traits_2    traits;

    Offset_polygon_with_holes_2                                 offset_conv;

    std::cout << "Using the convolution method ... ";
    offset_conv = offset_polygon_2 (pgn, r, traits);
    std::cout << "Done." << std::endl;

    // Define auxiliary polygon-decomposition objects.
    CGAL::Small_side_angle_bisector_decomposition_2<Rat_kernel> ssab_decomp;
    CGAL::Optimal_convex_decomposition_2<Rat_kernel>            opt_decomp;
    CGAL::Hertel_Mehlhorn_convex_decomposition_2<Rat_kernel>    hm_approx_decomp;
    CGAL::Greene_convex_decomposition_2<Rat_kernel>             greene_decomp;
    Offset_polygon_with_holes_2                                 offset_decomp;

    if (use_ssab)
    {
      std::cout << "Using the small-side angle-bisector decomposition ... ";
      offset_decomp = offset_polygon_2 (pgn, r, ssab_decomp, traits);
      if (are_equal (offset_conv, offset_decomp))
      {
        std::cout << "OK." << std::endl;
      }
      else
      {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (use_opt)
    {
      std::cout << "Using the optimal convex decomposition ... ";
      offset_decomp = offset_polygon_2 (pgn, r, opt_decomp, traits);
      if (are_equal (offset_conv, offset_decomp))
      {
        std::cout << "OK." << std::endl;
      }
      else
      {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (use_hm)
    {
      std::cout << "Using the Hertel--Mehlhorn decomposition ... ";
      offset_decomp = offset_polygon_2 (pgn, r, hm_approx_decomp, traits);
      if (are_equal (offset_conv, offset_decomp))
      {
        std::cout << "OK." << std::endl;
      }
      else
      {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (use_greene)
    {
      std::cout << "Using the Greene decomposition ... ";
      offset_decomp = offset_polygon_2 (pgn, r, greene_decomp, traits);
      if (are_equal (offset_conv, offset_decomp))
      {
        std::cout << "OK." << std::endl;
      }
      else
      {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (i+2 < argc && argv[i+2][0] == '-')
      i += 3;
    else
      i += 2;

  }
  return (0);
}

#endif
