//! \file examples/Arrangement_on_surface_2/dual_lines.cpp
// Checking whether there are three collinear points in a given input set
// using the arrangement of the dual lines.

#include <algorithm>

#include <CGAL/basic.h>
#include <CGAL/Arr_curve_data_traits_2.h>

#include "arr_linear.h"
#include "read_objects.h"

using Data_traits = CGAL::Arr_curve_data_traits_2<Traits, size_t>;
using Data_x_monotone_curve_2 = Data_traits::X_monotone_curve_2;
using Data_arrangement = CGAL::Arrangement_2<Data_traits>;

int main(int argc, char* argv[]) {
  // Get the name of the input file from the command line, or use the default
  // points.dat file if no command-line parameters are given.
  const char* filename = (argc > 1) ? argv[1] : "coll_points.dat";

  std::vector<Point> points;
  read_objects<Point>(filename, std::back_inserter(points));
  std::vector<Data_x_monotone_curve_2> dual_lines(points.size());
  size_t k{0};
  std::transform(points.begin(), points.end(), dual_lines.begin(),
                 [&](const Point& p) {
                   Line dual_line(p.x(), -1, -(p.y()));
                   return Data_x_monotone_curve_2(dual_line, k++);
                 });

  // Construct the dual arrangement by aggregately inserting the lines.
  Data_arrangement arr;
  insert(arr, dual_lines.begin(), dual_lines.end());

  // Look for vertices whose degree is greater than 4.
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->degree() > 4) {
      // There should be vit->degree()/2 lines intersecting at the current
      // vertex. We print their primal points and their indices.
      auto circ = vit->incident_halfedges();
      for (size_t d = 0; d < vit->degree() / 2; ++d) {
        k = circ->curve().data();     // The index of the primal point.
        std::cout << "Point no. " << k+1 << ": (" << points[k] << "), ";
        ++circ;
      }
      std::cout << "are collinear.\n";
    }
  }
  return 0;
}
