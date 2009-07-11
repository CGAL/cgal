//! \file examples/Arrangement_on_surface_2/dual_lines.cpp
// Checking whether there are three collinear points in a given input set
// using the arrangement of the dual lines.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <cstdlib>

typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>             Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Line_2                              Line_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main (int argc, char *argv[])
{
  // Get the name of the input file from the command line, or use the default
  // points.dat file if no command-line parameters are given.
  const char * filename = (argc > 1) ? argv[1] : "points.dat";

  // Open the input file.
  std::ifstream     in_file (filename);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << "!" << std::endl;
    return 1;
  }

  // Read the points from the file, and consturct their dual lines.
  // The input file format should be (all coordinate values are integers):
  // <n>                                 // number of point.
  // <x_1> <y_1>                         // point #1.
  // <x_2> <y_2>                         // point #2.
  //   :      :       :      :
  // <x_n> <y_n>                         // point #n.
  std::vector<Point_2>           points;
  std::list<X_monotone_curve_2>  dual_lines;

  unsigned int n;
  in_file >> n;
  points.resize(n);
  unsigned int k;
  for (k = 0; k < n; ++k) {
    int px, py;
    in_file >> px >> py;
    points[k] = Point_2 (px, py);

    // The line dual to the point (p_x, p_y) is y = p_x*x - p_y,
    // or: p_x*x - y - p_y = 0:
    dual_lines.push_back (Line_2 (Number_type(px),
                                  Number_type(-1), 
                                  Number_type(-py)));
  }
  in_file.close();

  // Construct the dual arrangement by aggragately inserting the lines.
  Arrangement_2      arr;

  insert (arr, dual_lines.begin(), dual_lines.end());

  std::cout << "The dual arrangement size:" << std::endl
            << "V = " << arr.number_of_vertices()
            << " (+ " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces()
            << " (" << arr.number_of_unbounded_faces()
            << " unbounded)" << std::endl;

  // Look for a vertex whose degree is greater than 4.
  Arrangement_2::Vertex_const_iterator   vit;
  bool                                   found_collinear = false;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->degree() > 4) {
      found_collinear = true;
      break;
    }
  }
  if (found_collinear)
    std::cout << "Found at least three collinear points in the input set."
              << std::endl;
  else
    std::cout << "No three collinear points are found in the input set."
              << std::endl;

  // Pick two points from the input set, compute their midpoint and insert
  // its dual line into the arrangement.
  Kernel             ker;
  const int          k1 = std::rand() % n, k2 = (k1 + 1) % n;
  Point_2            p_mid = ker.construct_midpoint_2_object() (points[k1],
                                                                points[k2]);
  X_monotone_curve_2 dual_p_mid = Line_2 (Number_type(p_mid.x()),
                                          Number_type(-1), 
                                          Number_type(-p_mid.y()));

  insert (arr, dual_p_mid);

  // Make sure that we now have three collinear points.
  found_collinear = false;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->degree() > 4) {
      found_collinear = true;
      break;
    }
  }
  CGAL_assertion (found_collinear);
  return (0);
}
