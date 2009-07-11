//! \file examples/Arrangement_on_surface_2/dual_lines.cpp
// Checking whether there are three collinear points in a given input set
// using the arrangement of the dual lines.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Cartesian<Number_type>                     Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>                Linear_traits_2;
typedef Linear_traits_2::Point_2                         Point_2;
typedef Linear_traits_2::Line_2                          Line_2;
typedef CGAL::Arr_curve_data_traits_2<Linear_traits_2,
                                      unsigned int>      Traits_2; 
typedef Traits_2::X_monotone_curve_2                     X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                    Arrangement_2;

int main (int argc, char *argv[])
{
  // Get the name of the input file from the command line, or use the default
  // points.dat file if no command-line parameters are given.
  const char * filename = (argc > 1) ? argv[1] : "coll_points.dat";

  // Open the input file.
  std::ifstream     in_file (filename);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << " ..." << std::endl;
    return (1);
  }

  // Read the points from the file, and consturct their dual lines.
  std::vector<Point_2>           points;
  std::list<X_monotone_curve_2>  dual_lines;

  unsigned int n;
  in_file >> n;
  points.resize (n);
  unsigned int k;
  for (k = 0; k < n; ++k) {
    int px, py;
    in_file >> px >> py;
    points[k] = Point_2 (px, py);

    // The line dual to the point (p_x, p_y) is y = p_x*x - p_y,
    // or: p_x*x - y - p_y = 0:
    Line_2 dual_line = Line_2(Number_type(px),
                              Number_type(-1),
                              Number_type(-py));

    // Generate the x-monotone curve based on the line and the point index.
    dual_lines.push_back (X_monotone_curve_2 (dual_line, k));
  }
  in_file.close();

  // Construct the dual arrangement by aggragately inserting the lines.
  Arrangement_2      arr;

  insert (arr, dual_lines.begin(), dual_lines.end());

  // Look for vertices whose degree is greater than 4.
  Arrangement_2::Vertex_const_iterator                    vit;
  Arrangement_2::Halfedge_around_vertex_const_circulator  circ;
  unsigned int                                            d;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->degree() > 4) {
      // There should be vit->degree()/2 lines intersecting at the current
      // vertex. We print their primal points and their indices.
      circ = vit->incident_halfedges();
      for (d = 0; d < vit->degree() / 2; d++) {
        k = circ->curve().data();     // The index of the primal point.
        std::cout << "Point no. " << k+1 << ": (" << points[k] << "), ";
        ++circ;
      }
      std::cout << "are collinear." << std::endl;
    }
  }
  return 0;
}
