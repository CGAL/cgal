//! \file examples/Arrangement_2/ex_Bezier_curves.C
// Constructing an arrangement of Bezier curves.
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

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             NT;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel,
                                        Alg_kernel,
                                        Nt_traits>      Traits_2;
typedef Traits_2::Curve_2                               Bezier_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

bool read_Bezier_curves (const char* filename,
                         std::list<Bezier_curve_2>& curves)
{
  // Clear the output list.
  curves.clear();

  // Open the input file.
  std::ifstream   in_file (filename);

  if (! in_file.is_open())
    return (false);
  
  // Read the number of curves.
  int             n_curves = -1;

  in_file >> n_curves;
  if (n_curves <= 0)
    return (false);

  // Read the curves one by one. Each curve is given by a number of control
  // points, followed by a list of the control points (with rational
  // coordinates).
  int             n_points;
  Rational        x, y;
  int             i, j;
  
  for (i = 0; i < n_curves; i++)
  {
    // Read the number of control points.
    in_file >> n_points;
    if (n_points <= 0)
      return (false);

    // Read the control points.
    std::vector<Rat_point_2>    ctrl_pts (n_points);

    for (j = 0; j < n_points; j++)
    {
      in_file >> x >> y;
      ctrl_pts[j] = Rat_point_2 (x, y);
    }

    // Construct the Bezier curve.
    curves.push_back (Bezier_curve_2 (ctrl_pts.begin(), ctrl_pts.end()));
  }

  return (true);
}

int main (int argc, char **argv)
{
  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.
  char   *filename = "Bezier.dat";

  if (argc > 1)
    filename = argv[1];

  // Read the curves from the input file.
  std::list<Bezier_curve_2>  curves;

  if (! read_Bezier_curves (filename, curves))
  {
    std::cerr << "Failed to read data from " << filename << std::endl;
    return (1);
  }

  // Print the curves.
  std::list<Bezier_curve_2>::iterator    cit;
  unsigned int                           ind = 1;

  for (cit = curves.begin(); cit != curves.end(); ++cit, ind++)
    std::cout << "B_" << ind << " (t) = " << *cit << std::endl;

  // Construct the arrangement.
  Arrangement_2                     arr;

  insert_curves (arr, curves.begin(), curves.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  return (0);
}

#endif

