//! \file examples/Arrangement_on_surface_2/Bezier_curves.cpp
// Constructing an arrangement of Bezier curves.

#include <CGAL/config.h>

#ifdef CGAL_USE_CORE

#include "arr_Bezier.h"
#include "arr_print.h"
#include "read_objects.h"

int main(int argc, char* argv[]) {
  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.
  const char* filename = (argc > 1) ? argv[1] : "Bezier.dat";

  // Read the Bezier curves.
  std::list<Bezier_curve>  curves;
  read_objects<Bezier_curve>(filename, std::back_inserter(curves));

  // Construct the arrangement.
  Arrangement arr;
  CGAL::insert(arr, curves.begin(), curves.end());
  print_arrangement_size(arr);

  return 0;
}


#else

#include <iostream>

int main ()
{
  std::cout << "Sorry, this example needs GMP and CORE\n";
  return 0;
}

#endif
