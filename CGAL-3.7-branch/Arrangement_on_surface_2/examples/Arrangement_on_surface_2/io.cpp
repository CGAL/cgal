//! \file examples/Arrangement_on_surface_2/io.cpp
// Using the arrangement I/O operators.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <fstream>

#include "point_location_utils.h"

typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
  // Construct the arrangement.
  Arrangement_2    arr;

  construct_segments_arr (arr);

  std::cout << "Writing an arrangement of size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Write the arrangement to a file.
  std::ofstream    out_file ("arr_ex_io.dat");

  out_file << arr;
  out_file.close();

  // Read the arrangement from the file.
  Arrangement_2    arr2;
  std::ifstream    in_file ("arr_ex_io.dat");

  in_file >> arr2;
  in_file.close();

  std::cout << "Read an arrangement of size:" << std::endl
            << "   V = " << arr2.number_of_vertices()
            << ",  E = " << arr2.number_of_edges()
            << ",  F = " << arr2.number_of_faces() << std::endl;

  return (0);
}
