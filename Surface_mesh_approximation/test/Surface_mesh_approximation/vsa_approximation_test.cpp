#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_approximation/approximate_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

/**
 * This file tests the free function CGAL::Surface_mesh_approximation::approximate_mesh.
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
  std::vector<Kernel::Vector_3> proxies;
  std::vector<Kernel::Point_3> points;
  std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles;

  CGAL::Surface_mesh_approximation::approximate_mesh(mesh,
    CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::INCREMENTAL).
      max_number_of_proxies(6).
      number_of_iterations(30).
      number_of_relaxations(5).
      subdivision_ratio(0.5).
      face_proxy_map(fpxmap).
      proxies(std::back_inserter(proxies)).
      anchors(std::back_inserter(points)).
      triangles(std::back_inserter(triangles)));

  std::cout << "#fpxmap " << fidxmap.size() << std::endl;
  std::cout << "#proxies " << proxies.size() << std::endl;
  std::cout << "#vertices " << points.size() << std::endl;
  std::cout << "#triangles " << triangles.size() << std::endl;

  return EXIT_SUCCESS;
}
