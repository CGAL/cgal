//! \file examples/Minkowski_sum_2/exact_offset.cpp
// Computing the exact offset of a polygon.
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
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/offset_polygon_2.h>
#include <CGAL/Timer.h>
#include <iostream>

typedef CGAL::CORE_algebraic_number_traits     Nt_traits;
typedef Nt_traits::Rational                    Rational;
typedef Nt_traits::Algebraic                   Algebraic;

struct Rat_kernel : public CGAL::Cartesian<Rational> {};
struct Alg_kernel : public CGAL::Cartesian<Algebraic> {};
struct Conic_traits_2 : public CGAL::Arr_conic_traits_2<Rat_kernel,
                                                        Alg_kernel,
                                                        Nt_traits> {};

typedef CGAL::Polygon_2<Rat_kernel>            Polygon_2;

typedef CGAL::Gps_traits_2<Conic_traits_2>     Gps_traits_2;
typedef Gps_traits_2::Polygon_2                Offset_polygon_2;
typedef Gps_traits_2::Polygon_with_holes_2     Offset_polygon_with_holes_2;

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

  // Compute the offset polygon.
  Conic_traits_2               traits;
  const Rational               radius = 5;
  Offset_polygon_with_holes_2  offset;
  CGAL::Timer                  timer;

  timer.start();
  offset = offset_polygon_2 (P, radius, traits);
  timer.stop();

  std::cout << "The offset polygon has "
            << offset.outer_boundary().size() << " vertices, "
            << offset.number_of_holes() << " holes." << std::endl;
  std::cout << "Offset computation took "
            << timer.time() << " seconds." << std::endl;
  return (0);
}

#endif
