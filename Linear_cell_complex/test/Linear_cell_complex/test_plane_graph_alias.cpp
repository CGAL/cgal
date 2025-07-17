#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <cassert>
#include <sstream>
#include <cstdlib>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2> LCC;

#ifndef CGAL_NO_DEPRECATED_CODE

int main()
{
  LCC lcc1, lcc2;

  // Planar graph with 3 vertices and 3 edges (triangle)
  std::stringstream input;
  input << "3 3\n0.0 0.0\n1.0 0.0\n0.0 1.0\n0 1\n1 2\n2 0\n";

  // Test new function
  auto d1 = CGAL::read_plane_graph_in_lcc(lcc1, input);
  assert(d1 != LCC::null_descriptor);

  // Rewind input stream for second test
  input.clear();
  input.seekg(0, std::ios::beg);

  // Test deprecated function
  auto d2 = CGAL::import_from_plane_graph(lcc2, input);
  assert(d2 != LCC::null_descriptor);

  return EXIT_SUCCESS;
}

#endif // CGAL_NO_DEPRECATED_CODE