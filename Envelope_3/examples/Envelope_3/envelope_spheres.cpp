//! \file examples/Envelope_3/envelope_spheres.cpp
// Constructing the lower envelope of a set of spheres read from a file.

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
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Env_sphere_traits_3.h>
#include <CGAL/envelope_3.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <list>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_3                           Rat_point_3;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;

typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                      Conic_traits_2;

typedef CGAL::Env_sphere_traits_3<Conic_traits_2>     Traits_3;
typedef Traits_3::Surface_3                           Sphere_3;
typedef CGAL::Envelope_diagram_2<Traits_3>            Envelope_diagram_2;

int main(int argc, char **argv)
{
  // Get the name of the input file from the command line, or use the default
  // fan_grids.dat file if no command-line parameters are given.
  const char * filename = (argc > 1) ? argv[1] : "spheres.dat";

  // Open the input file.
  std::ifstream     in_file(filename);

  if (! in_file.is_open())
  {
    std::cerr << "Failed to open " << filename << " ..." << std::endl;
    return 1;
  }

  // Read the spheres from the file.
  // The input file format should be (all coordinate values are integers):
  // <n>                           // number of spheres.
  // <x_1> <y_1> <x_1> <R_1>       // center and squared radious of sphere #1.
  // <x_2> <y_2> <x_2> <R_2>       // center and squared radious of sphere #2.
  //   :     :     :     :
  // <x_n> <y_n> <x_n> <R_n>       // center and squared radious of sphere #n.
  int                   n = 0;
  std::list<Sphere_3>   spheres;
  int                   x = 0, y = 0, z = 0, sqr_r = 0;
  int                   i;

  in_file >> n;
  for (i = 0; i < n; ++i)
  {
    in_file >> x >> y >> z >> sqr_r;
    spheres.push_back(Sphere_3(Rat_point_3(x, y, z), Rational(sqr_r)));
  }
  in_file.close();

  // Compute the lower envelope.
  Envelope_diagram_2    min_diag;
  CGAL::Timer           timer;

  std::cout << "Constructing the lower envelope of "
            << n << " spheres." << std::endl;

  timer.start();
  CGAL::lower_envelope_3(spheres.begin(), spheres.end(), min_diag);
  timer.stop();

  // Print the dimensions of the minimization diagram.
  std::cout << "V = " << min_diag.number_of_vertices()
            << ",  E = " << min_diag.number_of_edges()
            << ",  F = " << min_diag.number_of_faces() << std::endl;

  std::cout << "Construction took " << timer.time()
            << " seconds." << std::endl;

  return 0;
}

#endif
