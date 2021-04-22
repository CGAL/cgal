// Constructing a arrangement of Bezier polycurves.

#include <CGAL/config.h>

#ifndef CGAL_USE_CORE

#include <iostream>

int main() {
  std::cout << "Sorry, this example needs CORE ...\n";
  return 0;
}

#else

#include <CGAL/basic.h>
#include <CGAL/Arr_polycurve_traits_2.h>

#include "arr_Bezier.h"
#include "arr_print.h"

typedef CGAL::Arr_polycurve_traits_2<Traits>            Polycurve_bezier_traits;
typedef Polycurve_bezier_traits::Point_2                Point;
typedef Polycurve_bezier_traits::X_monotone_curve_2     X_mono_polycurve;
typedef CGAL::Arrangement_2<Polycurve_bezier_traits>    Arrangement_2;

typedef boost::variant<Point, Bezier_x_monotone_curve>  Make_x_monotone_result;

int main() {
  Polycurve_bezier_traits pc_traits;
  Traits bezier_traits;

  auto ctr_xpolycurve = pc_traits.construct_x_monotone_curve_2_object();

  std::vector<Bezier_x_monotone_curve> x_bezier_curves;

  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.
  const char* filename = "Bezier_polycurve.dat";

  // Open the input file.
  std::ifstream in_file(filename);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return 1;
  }

  // Read the curves from the input file.
  unsigned int n_curves;
  std::list<Bezier_x_monotone_curve> x_curves;
  Bezier_curve B;

  in_file >> n_curves;
  for (size_t k = 0; k < n_curves; ++k) {
    // Read the current curve (specified by its control points).
    in_file >> B;
    // convert it into x-monotone bezier curve.
    std::vector<Make_x_monotone_result> obj_vector;
    bezier_traits.make_x_monotone_2_object()(B, std::back_inserter(obj_vector));
    auto* x_seg_p = boost::get<Bezier_x_monotone_curve>(&obj_vector[0]);
    CGAL_assertion(x_seg_p);
    x_curves.push_back(*x_seg_p);
  }

  X_mono_polycurve polycurve = ctr_xpolycurve(x_curves.begin(), x_curves.end());
  Arrangement_2 arr;
  insert(arr, polycurve);               // construct the arrangement
  print_arrangement_size(arr);          // print the arrangement size

  return 0;
}

#endif
