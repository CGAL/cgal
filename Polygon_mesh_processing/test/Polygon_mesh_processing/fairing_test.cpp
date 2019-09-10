#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polygon_mesh_processing/fair.h>

#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>
#include <iterator>

typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;

template <typename K>
void test_polyhedron(const char* filename, const K&, const bool save_output)
{
  typedef CGAL::Polyhedron_3<K> Polyhedron;

  //run test for a Polyhedron
  Polyhedron poly; // file should contain oriented polyhedron
  std::ifstream input(filename);

  if (!input || !(input >> poly))
  {
    std::cerr << "Error: cannot read Polyhedron : " << filename << "\n";
    assert(!poly.empty());
    assert(false);
    return;
  }

  //try to fair the mesh
  std::size_t nbv =
    std::distance(vertices(poly).first, vertices(poly).second);

  CGAL::Polygon_mesh_processing::fair(poly, vertices(poly));

  std::size_t nbv2 =
    std::distance(vertices(poly).first, vertices(poly).second);

  assert(nbv == nbv2);

  if (!save_output)
    return;

  std::ofstream faired_off("faired.off");
  faired_off << poly;
  faired_off.close();
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/elephant.off";
  const bool save_output = (argc > 2) ? true : false;

  test_polyhedron(filename, Epic(), save_output);
  test_polyhedron(filename, Epec(), save_output);

  std::cerr << "All done." << std::endl;

  return 0;
}
