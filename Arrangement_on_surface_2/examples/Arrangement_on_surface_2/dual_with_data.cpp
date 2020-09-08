//! \file examples/Arrangement_on_surface_2/dual_lines.cpp
// Checking whether there are three collinear points in a given input set
// using the arrangement of the dual lines.

#include <algorithm>

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include "read_objects.h"

typedef CGAL::Cartesian<CGAL::Exact_rational>            Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>                Linear_traits_2;
typedef Linear_traits_2::Point_2                         Point_2;
typedef Linear_traits_2::Line_2                          Line_2;
typedef CGAL::Arr_curve_data_traits_2<Linear_traits_2,
                                      unsigned int>      Traits_2;
typedef Traits_2::X_monotone_curve_2                     X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                    Arrangement_2;

int main(int argc, char *argv[])
{
  // Get the name of the input file from the command line, or use the default
  // points.dat file if no command-line parameters are given.
  const char* filename = (argc > 1) ? argv[1] : "coll_points.dat";

  std::vector<Point_2> points;

  read_objects<Point_2>(filename, std::back_inserter(points));
  std::vector<X_monotone_curve_2> dual_lines(points.size());
  size_t k{0};
  std::transform(points.begin(), points.end(), dual_lines.begin(),
                 [&](const Point_2& p) {
                   Line_2 dual_line(p.x(), CGAL::Exact_rational(-1), -(p.y()));
                   return X_monotone_curve_2(dual_line, k++);
                 });

  // Construct the dual arrangement by aggregately inserting the lines.
  Arrangement_2 arr;
  insert(arr, dual_lines.begin(), dual_lines.end());

  // Look for vertices whose degree is greater than 4.
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->degree() > 4) {
      // There should be vit->degree()/2 lines intersecting at the current
      // vertex. We print their primal points and their indices.
      auto circ = vit->incident_halfedges();
      for (auto d = 0; d < vit->degree() / 2; ++d) {
        k = circ->curve().data();     // The index of the primal point.
        std::cout << "Point no. " << k+1 << ": (" << points[k] << "), ";
        ++circ;
      }
      std::cout << "are collinear." << std::endl;
    }
  }
  return 0;
}
