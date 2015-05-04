#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polygon_mesh_processing/fair.h>

#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>
#include <iterator>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::face_descriptor   face_descriptor;

void test(const char* file_name)
{
  //run test for a Polyhedron
  std::ifstream input(file_name);
  Polyhedron poly; // file should contain oriented polyhedron
  if (!input || !(input >> poly) || poly.empty())
  {
    std::cerr << "Error: cannot read Polyhedron : " << file_name << "\n";
    CGAL_assertion(false);
  }

  //try to fair the mesh
  unsigned int nbv = std::distance(vertices(poly).first
                                 , vertices(poly).second);

  CGAL::Polygon_mesh_processing::fair(poly, vertices(poly));

  unsigned int nbv2 = std::distance(vertices(poly).first
                                  , vertices(poly).second);

  CGAL_assertion(nbv == nbv2);

  std::ofstream faired_off("faired.off");
  faired_off << poly;
  faired_off.close();
}

int main()
{
  test("data/elephant.off");

  std::cerr << "All done." << std::endl;
}
