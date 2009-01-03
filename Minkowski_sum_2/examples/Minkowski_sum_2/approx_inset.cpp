//! \file examples/Minkowski_sum_2/approx_inset.cpp
// Computing the approximated inset of a polygon.

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
  std::ifstream    in_file ("tight.dat");

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
  const Number_type                      radius = 1;
  const double                           err_bound = 0.00001;
  std::list<Offset_polygon_2>            inset_polygons;
  std::list<Offset_polygon_2>::iterator  iit;
  CGAL::Timer                            timer;

  timer.start();
  approximated_inset_2 (P, radius, err_bound,
                        std::back_inserter (inset_polygons));
  timer.stop();

  std::cout << "The inset comprises " 
            << inset_polygons.size() << " polygon(s)." << std::endl;
  for (iit = inset_polygons.begin(); iit != inset_polygons.end(); ++iit)
  {
      std::cout << "    Polygon with "
                << iit->size() << " vertices." << std::endl;
  }
  std::cout << "Inset computation took "
            << timer.time() << " seconds." << std::endl;
  return (0);
}
