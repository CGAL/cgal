//! \file examples/Arrangement_2/io_unbounded.cpp
// Using the I/O operators with an arrangement of unbounded curves.

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <list>
#include <fstream>

typedef CGAL::Cartesian<CGAL::Exact_rational>         Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>             Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Segment_2                           Segment_2;
typedef Traits_2::Ray_2                               Ray_2;
typedef Traits_2::Line_2                              Line_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
  // Construct an arrangement of five linear objects.
  Arrangement_2                  arr;
  std::list<X_monotone_curve_2>  curves;

  curves.push_back (Line_2 (Point_2 (0, 0), Point_2 (2, 1)));
  curves.push_back (Line_2 (Point_2 (0, 0), Point_2 (2, -1)));
  curves.push_back (Line_2 (Point_2 (-1, 0), Point_2 (-1, 1)));
  curves.push_back (Ray_2 (Point_2 (2, 3), Point_2 (2, 4)));
  curves.push_back (Segment_2 (Point_2 (0, 1), Point_2 (0, 2)));

  insert (arr, curves.begin(), curves.end());

  // Print out the size of the resulting arrangement.
  std::cout << "Writing an arrangement of size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << " (plus " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces()
            << " (" << arr.number_of_unbounded_faces() << " unbounded)"
            << std::endl << std::endl;

  // Write the arrangement to a file.
  std::ofstream    out_file ("arr_ex_io_unbounded.dat");

  out_file << arr;
  out_file.close();

  // Read the arrangement from the file.
  Arrangement_2    arr2;
  std::ifstream    in_file ("arr_ex_io_unbounded.dat");

  in_file >> arr2;
  in_file.close();

  std::cout << "Read an arrangement of size:" << std::endl
            << "   V = " << arr2.number_of_vertices()
            << " (plus " << arr2.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr2.number_of_edges()
            << ",  F = " << arr2.number_of_faces()
            << " (" << arr2.number_of_unbounded_faces() << " unbounded)"
            << std::endl << std::endl;

  return (0);
}
