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
  std::vector<std::vector<std::size_t> > triangles;
  std::vector<Kernel::Point_3> points;
  std::list<Polyhedron::Vertex_handle> anchors;
  std::vector<CGAL::VSA::Plane_proxy<Kernel> > proxies;
  std::map<Polyhedron::Facet_handle, std::size_t> fidxmap;
  boost::associative_property_map<std::map<Polyhedron::Facet_handle, std::size_t> > fpxmap(fidxmap);

  CGAL::VSA::mesh_approximation(mesh,
    std::back_inserter(points),
    std::back_inserter(triangles),
    CGAL::VSA::parameters::max_nb_proxies(6).
      nb_of_iterations(30).
      nb_of_relaxations(5).
      mesh_chord_error(0.5).
      facet_proxy_map(fpxmap).
      anchor_vertices(std::back_inserter(anchors)).
      proxies(std::back_inserter(proxies)).
      output_mesh(&out_mesh));

  std::cout << "#triangles " << triangles.size() << std::endl;
  std::cout << "#vertices " << points.size() << std::endl;
  std::cout << "#anchor_vertices " << anchors.size() << std::endl;
  std::cout << "#proxies " << proxies.size() << std::endl;
  std::cout << "#fpxmap " << fidxmap.size() << std::endl;

  return EXIT_SUCCESS;
}
