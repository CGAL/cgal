#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

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
  typedef CGAL::Surface_mesh<typename K::Point_3> Polyhedron;
  typedef typename Polyhedron::Vertex_index Vertex_index;

  //run test for a Polyhedron
  Polyhedron poly; // file should contain oriented polyhedron
  std::ifstream input(filename);

  if (!input || !(input >> poly))
  {
    std::cerr << "Error: cannot read Polyhedron : " << filename << "\n";
    assert(false);
    return;
  }
  assert(!poly.is_empty());

  //try to fair the mesh
  std::size_t nbv =
    std::distance(vertices(poly).first, vertices(poly).second);
  std::vector<typename Polyhedron::Vertex_index> sel_vert(10);
  sel_vert[0] = Vertex_index(142);
  sel_vert[1] = Vertex_index(237);
  sel_vert[2] = Vertex_index(1499);
  sel_vert[3] = Vertex_index(1498);
  sel_vert[4] = Vertex_index(74);
  sel_vert[5] = Vertex_index(141);
  sel_vert[6] = Vertex_index(2054);
  sel_vert[7] = Vertex_index(2053);
  sel_vert[8] = Vertex_index(140);
  sel_vert[9] = Vertex_index(2052);
  CGAL::Polygon_mesh_processing::fair(poly, sel_vert);

  std::size_t nbv2 =
    std::distance(vertices(poly).first, vertices(poly).second);

  assert(nbv == nbv2);

  if (!save_output)
    return;

  std::ofstream faired_off("faired.off");
  faired_off << poly;
  faired_off.close();
}

int main()
{
  const char* filename = "data/elephant.off";
    test_polyhedron(filename, Epic(), false);
    test_polyhedron(filename, Epec(), false);

  std::cerr << "All done." << std::endl;

  return 0;
}
