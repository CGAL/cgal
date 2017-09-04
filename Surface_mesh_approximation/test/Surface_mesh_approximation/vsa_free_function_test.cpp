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
  std::vector<CGAL::PlaneProxy<Kernel> > proxies;
  std::map<Polyhedron::Facet_handle, std::size_t> fidxmap;
  boost::associative_property_map<std::map<Polyhedron::Facet_handle, std::size_t> > fpxmap(fidxmap);
  CGAL::vsa_mesh_approximation(mesh, out_mesh,
    CGAL::VSA::parameters::seeding_by_number(6).
      iterations(30).
      inner_iterations(5).
      chord_subdivide(0.5).
      facet_proxy_map(fpxmap).
      anchor_vertex(std::back_inserter(anchor_vtx)).
      anchor_point(std::back_inserter(anchor_pos)).
      indexed_triangles(std::back_inserter(tris)).
      proxies(std::back_inserter(proxies)));

  std::cout << "#tris " << tris.size() << std::endl;
  std::cout << "#anchor_pos " << anchor_pos.size() << std::endl;
  std::cout << "#anchor_vtx " << anchor_vtx.size() << std::endl;
  std::cout << "#proxies " << proxies.size() << std::endl;
  std::cout << "#fpxmap " << fidxmap.size() << std::endl;

  return EXIT_SUCCESS;
}
