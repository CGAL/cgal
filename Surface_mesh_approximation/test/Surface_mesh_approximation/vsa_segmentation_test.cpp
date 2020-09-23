#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Mesh::Property_map<face_descriptor, std::size_t> Face_proxy_pmap;

/**
 * This file tests the free function CGAL::Surface_mesh_approximation::approximate_triangle_mesh.
 */
int main()
{
  Mesh mesh;
  std::ifstream file("data/sphere.off");
  if (!file || !(file >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  Face_proxy_pmap fpxmap =
    mesh.add_property_map<face_descriptor, std::size_t>("f:proxy_id", 0).first;
  std::vector<Kernel::Vector_3> proxies;

  // free function interface with named parameters
  CGAL::Surface_mesh_approximation::approximate_triangle_mesh(mesh,
    CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::HIERARCHICAL). // hierarchical seeding
    max_number_of_proxies(200). // both maximum number of proxies stop criterion,
    min_error_drop(0.05). // and minimum error drop stop criterion are specified
    number_of_iterations(30). // number of clustering iterations after seeding
    number_of_relaxations(5). // number of relaxations in seeding
    face_proxy_map(fpxmap). // output indexed face set
    proxies(std::back_inserter(proxies))); // number of iterations after seeding

  return EXIT_SUCCESS;
}
