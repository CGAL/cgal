//! \file examples/Minkowski_sum_2/approx_offset.cpp
// Computing the approximated offset of a polygon.

#include "ms_rational_nt.h"
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Cartesian.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/offset_polygon_2.h>
#include <CGAL/Timer.h>
#include <iostream>

typedef CGAL::Lazy_exact_nt<Number_type>           Lazy_exact_nt;

struct Kernel : public CGAL::Cartesian<Lazy_exact_nt> {};
typedef CGAL::Polygon_2<Kernel>                    Polygon_2;

typedef CGAL::Gps_circle_segment_traits_2<Kernel>  Gps_traits_2;
typedef Gps_traits_2::Polygon_2                    Offset_polygon_2;
typedef Gps_traits_2::Polygon_with_holes_2         Offset_polygon_with_holes_2;

int main ()
{
  // Open the input file.
  std::ifstream    in_file ("spiked.dat");

  if (! in_file.is_open())
  {
    std::cerr << "Failed to open the input file." << std::endl;
    return (1);
  }

  // Read the input polygon.
  Polygon_2        P;

  in_file >> P;
  in_file.close();

  std::cout << "Read an input polygon with "
            << P.size() << " vertices." << std::endl;

  // Approximate the offset polygon.
  const Number_type            radius = 5;
  const double                 err_bound = 0.00001;
  Offset_polygon_with_holes_2  offset;
  CGAL::Timer                  timer;

  timer.start();
  offset = approximated_offset_2 (P, radius, err_bound);
  timer.stop();

  std::cout << "The offset polygon has "
            << offset.outer_boundary().size() << " vertices, "
            << offset.number_of_holes() << " holes." << std::endl;
  std::cout << "Offset computation took "
            << timer.time() << " seconds." << std::endl;
  return (0);
}
