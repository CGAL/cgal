//! \file examples/Arrangement_on_surface_2/predefined_kernel_non_intersecting.cpp
// Constructing an arrangement of non-intersecting line segments using the
// predefined kernel with exact predicates.

#include <list>
#include <fstream>
#include <chrono>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_non_caching_segment_basic_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Number_type = Kernel::FT;
using Traits_2 = CGAL::Arr_non_caching_segment_basic_traits_2<Kernel>;
using Point_2 = Traits_2::Point_2;
using Segment_2 = Traits_2::X_monotone_curve_2;
using Arrangement_2 = CGAL::Arrangement_2<Traits_2>;

int main(int argc, char* argv[]) {
  // Get the name of the input file from the command line, or use the default
  // Europe.dat file if no command-line parameters are given.
  const char* filename = (argc > 1) ? argv[1] : "Europe.dat";

  // Open the input file.
  std::ifstream in_file(filename);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << " ...\n";
    return 1;
  }

  // Read the segments from the file.
  // The input file format should be (all coordinate values are double
  // precision floating-point numbers):
  // <n>                                 // number of segments.
  // <sx_1> <sy_1>  <tx_1> <ty_1>        // source and target of segment #1.
  // <sx_2> <sy_2>  <tx_2> <ty_2>        // source and target of segment #2.
  //   :      :       :      :
  // <sx_n> <sy_n>  <tx_n> <ty_n>        // source and target of segment #n.
  std::list<Segment_2>  segments;

  std::size_t n;
  in_file >> n;
  for (std::size_t i = 0; i < n; ++i) {
    double sx, sy, tx, ty;
    in_file >> sx >> sy >> tx >> ty;
    segments.push_back(Segment_2(Point_2 (Number_type(sx), Number_type(sy)),
                                 Point_2 (Number_type(tx), Number_type(ty))));
  }
  in_file.close();

  // Construct the arrangement by aggregately inserting all segments.
  Arrangement_2 arr;

  std::cout << "Performing aggregated insertion of " << n << " segments.\n";

  auto start = std::chrono::system_clock::now();
  insert_non_intersecting_curves (arr, segments.begin(), segments.end());
  std::chrono::duration<double> secs = std::chrono::system_clock::now() - start;

  // Print the arrangement dimensions.
  std::cout << "V = " << arr.number_of_vertices()
            << ", E = " << arr.number_of_edges()
            << ", F = " << arr.number_of_faces() << std::endl;

  std::cout << "Construction took " << secs.count() << " seconds.\n";

  return 0;
}
