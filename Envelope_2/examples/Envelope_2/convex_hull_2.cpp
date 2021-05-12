//! \file examples/Envelope_2/ex_convex_hull.cpp
// Compute the convex hull of set of points using the lower envelope and upper
// envelopes of their dual line.

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>
#include <vector>

typedef CGAL::Exact_rational                             Number_type;
typedef CGAL::Cartesian<Number_type>                     Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>                Linear_traits_2;
typedef Linear_traits_2::Point_2                         Point_2;
typedef Linear_traits_2::Line_2                          Line_2;
typedef CGAL::Arr_curve_data_traits_2<Linear_traits_2,
                                      unsigned int>      Traits_2;
typedef Traits_2::X_monotone_curve_2                     Dual_line_2;
typedef CGAL::Envelope_diagram_1<Traits_2>               Diagram_1;

int main (int argc, char* argv[])
{
  // Read the points from the input file.
  const char* filename = (argc > 1) ? argv[1] : "ch_points.dat";
  std::ifstream in_file(filename);
  if (!in_file.is_open()) {
    std::cerr << "Failed to open " << filename << " ..." << std::endl;
    return -1;
  }

  // Read the points from the file, and construct their dual lines.
  std::list<Dual_line_2>  dual_lines;

  unsigned int            n;
  in_file >> n;
  std::vector<Point_2>    points;
  points.resize(n);

  for (unsigned int k = 0; k < n; ++k) {
    int px, py;
    in_file >> px >> py;
    points[k] = Point_2 (px, py);

    // The line dual to the point (p_x, p_y) is y = p_x*x - p_y,
    // or: p_x*x - y - p_y = 0:
    Line_2 line = Line_2 (Number_type(px), Number_type(-1), Number_type(-py));

    // Generate the x-monotone curve based on the line and the point index.
    dual_lines.push_back (Dual_line_2 (line, k));
  }
  in_file.close();

  // Compute the lower envelope of dual lines, which corresponds to the upper
  // part of the convex hull, and their upper envelope, which corresponds to
  // the lower part of the convex hull.
  Diagram_1              min_diag;
  Diagram_1              max_diag;
  lower_envelope_x_monotone_2(dual_lines.begin(), dual_lines.end(), min_diag);
  upper_envelope_x_monotone_2(dual_lines.begin(), dual_lines.end(), max_diag);

  // Output the points along the boundary convex hull in counterclockwise
  // order. We start by traversing the minimization diagram from left to
  // right, then the maximization diagram from right to left.
  std::cout << "The convex hull of " << points.size() << " input points:";
  Diagram_1::Edge_const_handle  e = min_diag.leftmost();
  while (e != min_diag.rightmost()) {
    std::cout << " (" << points[e->curve().data()] << ')';
    e = e->right()->right();
  }

  // Handle the degenerate case of a vertical convex hull edge:
  if (e->curve().data() != max_diag.leftmost()->curve().data())
    std::cout << " (" << points[e->curve().data()] << ')';

  e = max_diag.leftmost();
  while (e != max_diag.rightmost()) {
    std::cout << " (" << points[e->curve().data()] << ')';
    e = e->right()->right();
  }
  std::cout << std::endl;

  return 0;
}
