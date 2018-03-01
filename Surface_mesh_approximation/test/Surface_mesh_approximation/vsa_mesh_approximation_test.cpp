#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

/**
 * This file tests the free function CGAL::mesh_approximation.
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
  std::map<Polyhedron::Facet_handle, std::size_t> fidxmap;
  boost::associative_property_map<std::map<Polyhedron::Facet_handle, std::size_t> > fpxmap(fidxmap);
  std::vector<CGAL::Plane_proxy<Kernel> > proxies;
  std::vector<Kernel::Point_3> points;
  std::vector<std::vector<std::size_t> > triangles;

  CGAL::mesh_approximation(mesh,
    CGAL::Surface_mesh_approximation::parameters::seeding_method(CGAL::Incremental).
      max_nb_proxies(6).
      nb_of_iterations(30).
      nb_of_relaxations(5).
      mesh_chord_error(0.5).
      facet_proxy_map(fpxmap).
      proxies(std::back_inserter(proxies)).
      anchor_points(std::back_inserter(points)).
      indexed_triangles(std::back_inserter(triangles)));

  std::cout << "#fpxmap " << fidxmap.size() << std::endl;
  std::cout << "#proxies " << proxies.size() << std::endl;
  std::cout << "#vertices " << points.size() << std::endl;
  std::cout << "#triangles " << triangles.size() << std::endl;

  return EXIT_SUCCESS;
}
