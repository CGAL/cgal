#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

/**
 * This file tests the VSA free function API.
 */
int main()
{
  Polyhedron mesh;
  std::ifstream input("./data/cube_meshed.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  Polyhedron out_mesh;
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;
  std::list<Polyhedron::Vertex_handle> anchor_vtx;
  CGAL::vsa_mesh_approximation(mesh, out_mesh,
    CGAL::VSA::parameters::number_of_proxies(6).
      number_of_iterations(30).
      anchor_vertex(std::back_inserter(anchor_vtx)).
      anchor_point(std::back_inserter(anchor_pos)).
      indexed_triangles(std::back_inserter(tris)));

  return EXIT_SUCCESS;
}
