/*! \file bezier_traits_adapter.cpp
 * Using the traits adaptor to generate a traits class for Bezier polygons.
 */

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return (0);
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Timer.h>

#include <fstream>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;

typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Traits_2;
  
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef Traits_2::Curve_2                               Bezier_curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Gps_traits_2<Traits_2>                    Gps_traits_2;
typedef Gps_traits_2::General_polygon_2                 Polygon_2;
typedef Gps_traits_2::General_polygon_with_holes_2      Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                 Polygon_set;

/*! Read a general polygon with holes, formed by Bezier curves, from the
 * given input file.
 */
bool read_Bezier_polygon (const char* filename, Polygon_with_holes_2& P)
{
  // Open the input file.
  std::ifstream   in_file (filename);

  if (! in_file.is_open())
    return false;
  
  // Read the number of curves.
  unsigned int      n_curves;
  unsigned int      k;

  in_file >> n_curves;

  // Read the curves one by one, and construct the general polygon these
  // curve form (the outer boundary and the holes inside it).
  Traits_2                    traits;
  Traits_2::Make_x_monotone_2 mk_x_monotone = traits.make_x_monotone_2_object();
  bool                        first = true;
  Rat_point_2                 p_0;
  std::list<X_monotone_curve_2>  xcvs;
  Rat_kernel                  ker;
  Rat_kernel::Equal_2         equal = ker.equal_2_object();
  std::list<Polygon_2>        pgns;

  for (k = 0; k < n_curves; k++) {
    // Read the current curve and subdivide it into x-monotone subcurves.
    Bezier_curve_2                           B;
    std::list<CGAL::Object>                  x_objs;
    std::list<CGAL::Object>::const_iterator  xoit;
    X_monotone_curve_2                       xcv;

    in_file >> B;
    mk_x_monotone (B, std::back_inserter (x_objs));
    
    for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) {
      if (CGAL::assign (xcv, *xoit))
        xcvs.push_back (xcv);
    }
    
    // Check if the current curve closes a polygon, namely whether it target
    // point (the last control point) equals the source of the first curve in
    // the current chain.
    if (! first) {
      if (equal (p_0, B.control_point(B.number_of_control_points() - 1))) {
        // Push a new polygon into the polygon list. Make sure that the polygon
        // is counterclockwise oriented if it represents the outer boundary
        // and clockwise oriented if it represents a hole.
        Polygon_2          pgn (xcvs.begin(), xcvs.end());
        CGAL::Orientation  orient = pgn.orientation();
        
        if ((pgns.empty() && orient == CGAL::CLOCKWISE) ||
            (! pgns.empty() && orient == CGAL::COUNTERCLOCKWISE))
          pgn.reverse_orientation();
        
        pgns.push_back (pgn);
        xcvs.clear();
        first = true;
      }
    }
    else {
      // This is the first curve in the chain - store its source point.
      p_0 = B.control_point(0);
      first = false;
    }
  }

  if (! xcvs.empty())
    return false;

  // Construct the polygon with holes.
  std::list<Polygon_2>::iterator     pit = pgns.begin();
  
  ++pit;
  P = Polygon_with_holes_2 (pgns.front(), pit, pgns.end());

  return true;
}

// The main program.
int main (int argc, char* argv[])
{
  // Get the name of the input files from the command line, or use the default
  // char_g.dat and char_m.dat files if no command-line parameters are given.
  const char* filename1 = (argc > 1) ? argv[1] : "char_g.dat";
  const char* filename2 = (argc > 2) ? argv[2] : "char_m.dat";

  // Read the general polygons from the input files.
  CGAL::Timer           timer;
  Polygon_with_holes_2  P1, P2;

  timer.start();

  if (! read_Bezier_polygon (filename1, P1)) {
    std::cerr << "Failed to read " << filename1 << " ..." << std::endl;
    return 1;
  }

  if (! read_Bezier_polygon (filename2, P2)) {
    std::cerr << "Failed to read " << filename2 << " ..." << std::endl;
    return 1;
  }

  timer.stop();
  std::cout << "Constructed the input polygons in " << timer.time() 
            << " seconds." << std::endl << std::endl;
  
  // Compute the intersection of the two polygons.
  Polygon_set                     R;
  Polygon_set::const_iterator     rit;

  timer.reset();
  timer.start();
  CGAL::intersection (P1, P2, std::back_inserter(R));
  timer.stop();
  
  std::cout << "The intersection polygons are of sizes: {";
  for (rit = R.begin(); rit != R.end(); ++rit)
    std::cout << ' ' << rit->outer_boundary().size();
  std::cout << " }" << std::endl;
  std::cout << "The intersection computation took "
            << timer.time() << " seconds." << std::endl;

  return 0;
}

#endif
