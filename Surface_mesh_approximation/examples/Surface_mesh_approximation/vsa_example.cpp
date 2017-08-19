#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char *argv[])
{
  if (argc < 5)
    return EXIT_FAILURE;

  // create and read Polyhedron
  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  
  int init = std::atoi(argv[2]);
  if (init < 0 || init > 2)
    return EXIT_FAILURE;

  Polyhedron out_mesh;
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;
  std::list<Polyhedron::Vertex_handle> anchor_vtx;
  bool is_manifold = CGAL::vsa_mesh_approximation(mesh, out_mesh,
    CGAL::VSA::parameters::number_of_proxies(num_proxies).
      number_of_iterations(num_iterations).
      anchor_vertex(std::back_inserter(anchor_vtx)).
      anchor_point(std::back_inserter(anchor_pos)).
      indexed_triangles(std::back_inserter(tris)));

  std::cout << "#anchor_vtx " << anchor_vtx.size() << std::endl;
  std::cout << "#anchor_pos " << anchor_pos.size() << std::endl;
  std::cout << "#tris " << tris.size() << std::endl;

  return EXIT_SUCCESS;
}
