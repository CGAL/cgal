// Constructing a arrangement of Bezier polycurves.

#include <CGAL/config.h>

#ifndef CGAL_USE_CORE

#include <iostream>

int main()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_polycurve_traits_2.h>
#include "arr_print.h"

typedef CGAL::CORE_algebraic_number_traits             Nt_traits;
typedef Nt_traits::Rational                            NT;
typedef Nt_traits::Rational                            Rational;
typedef Nt_traits::Algebraic                           Algebraic;
typedef CGAL::Cartesian<Rational>                      Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                     Alg_kernel;
typedef Rat_kernel::Point_2                            Rat_point_2;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                       Bezier_traits;
typedef Bezier_traits::Curve_2                         Bezier_curve_2;
typedef Bezier_traits::X_monotone_curve_2              Bezier_x_curve_2;
typedef CGAL::Arr_polycurve_traits_2<Bezier_traits>    Polycurve_bezier_traits_2;
typedef Polycurve_bezier_traits_2::Point_2             Point_2;
typedef Polycurve_bezier_traits_2::X_monotone_curve_2  X_mono_polycurve;
typedef CGAL::Arrangement_2<Polycurve_bezier_traits_2> Arrangement_2;

typedef boost::variant<Point_2, Bezier_x_curve_2>      Make_x_monotone_result;

int main()
{
  Polycurve_bezier_traits_2 pc_traits;
  Bezier_traits bezier_traits;

  auto construct_x_mono_polycurve =
    pc_traits.construct_x_monotone_curve_2_object();

  std::vector<Bezier_x_curve_2> x_bezier_curves;
  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.
  const char* filename = "Bezier_polycurve.dat";

  // Open the input file.
  std::ifstream in_file (filename);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return 1;
  }

  // Read the curves from the input file.
  unsigned int n_curves;
  std::list<Bezier_x_curve_2> x_curves;
  Bezier_curve_2 B;

  in_file >> n_curves;
  unsigned int k;
  for (k = 0; k < n_curves; ++k) {
    // Read the current curve (specified by its control points).
    in_file >> B;
    //convert it into x-monotone bezier curve.
    std::vector<Make_x_monotone_result> obj_vector;
    bezier_traits.make_x_monotone_2_object()(B, std::back_inserter(obj_vector));
    auto* x_seg_p = boost::get<Bezier_x_curve_2>(&obj_vector[0]);
    CGAL_assertion(x_seg_p);
    x_curves.push_back(*x_seg_p);
  }

  X_mono_polycurve polycurve =
    construct_x_mono_polycurve(x_curves.begin(), x_curves.end());

  // Construct the arrangement.
  Arrangement_2 arr;
  insert(arr, polycurve);

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}
#endif
