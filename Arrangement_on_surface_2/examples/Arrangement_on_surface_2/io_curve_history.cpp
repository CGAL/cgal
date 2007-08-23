//! \file examples/Arrangement_on_surface_2/io_curve_history.cpp
// Using the arrangement-with-history I/O operators.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/IO/Arr_with_history_iostream.h>
#include <fstream>

typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Segment_2;
typedef CGAL::Arrangement_with_history_2<Traits_2>    Arr_with_hist_2;

int main ()
{
  Arr_with_hist_2   arr;

  // Insert six additional segments aggregately:
  Segment_2         segs[6];
  segs[0] = Segment_2 (Point_2 (2, 6), Point_2 (7, 1));
  segs[1] = Segment_2 (Point_2 (3, 2), Point_2 (3, 5));
  segs[2] = Segment_2 (Point_2 (2, 3), Point_2 (5, 3));
  segs[3] = Segment_2 (Point_2 (2, 6), Point_2 (7, 1));
  segs[4] = Segment_2 (Point_2 (0, 0), Point_2 (2, 6));
  segs[5] = Segment_2 (Point_2 (3, 4), Point_2 (6, 4));
  insert (arr, segs, segs + 6);

  std::cout << "Writing an arrangement of "
            << arr.number_of_curves() << " input segments:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Write the arrangement to a file.
  std::ofstream     out_file ("arr_ex_io_hist.dat");

  out_file << arr;
  out_file.close();

  // Read the arrangement from the file.
  Arr_with_hist_2   arr2;
  std::ifstream     in_file ("arr_ex_io_hist.dat");

  in_file >> arr2;
  in_file.close();

  std::cout << "Read an arrangement of "
            << arr2.number_of_curves() << " input segments:" << std::endl
            << "   V = " << arr2.number_of_vertices()
            << ",  E = " << arr2.number_of_edges()
            << ",  F = " << arr2.number_of_faces() << std::endl;

  return (0);
}
