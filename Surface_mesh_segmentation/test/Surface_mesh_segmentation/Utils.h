#include <fstream>
#include <iostream>
#include <CGAL/boost/graph/helpers.h>

template<class Polyhedron>
bool read_to_polyhedron(const char* file_name, Polyhedron& mesh)
{
  std::ifstream input(file_name);

  if (!input || !(input >> mesh))
  {
    std::cerr << "Failed to read mesh" << std::endl;
    return false;
  }

  if (CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Input mesh is invalid" << std::endl;
    return false;
  }

  return true;
}

void expect_or_fail(bool value)
{
  if (!value) exit(EXIT_FAILURE);
}
