#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

/**
 * This file tests the free function CGAL::Surface_mesh_approximation::approximate_triangle_mesh.
 */
int main()
{
  Mesh mesh;
  std::ifstream input("./data/cube.off");
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  const double target_edge_length = 0.05;
  const unsigned int nb_iter = 3;

  std::cout << "Start remeshing. ("
    << std::distance(faces(mesh).first, faces(mesh).second) << " faces)..." << std::endl;
  PMP::isotropic_remeshing(
    faces(mesh),
    target_edge_length,
    mesh,
    PMP::parameters::number_of_iterations(nb_iter));
  std::cout << "Remeshing done. ("
    << std::distance(faces(mesh).first, faces(mesh).second) << " faces)..." << std::endl;

  Mesh::Property_map<face_descriptor, std::size_t> fpxmap =
    mesh.add_property_map<face_descriptor, std::size_t>("f:proxy_id", 0).first;
  std::vector<Kernel::Vector_3> proxies;
  std::vector<Kernel::Point_3> points;
  std::vector<std::array<std::size_t, 3> > triangles;

  CGAL::Surface_mesh_approximation::approximate_triangle_mesh(mesh,
    CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::INCREMENTAL).
      max_number_of_proxies(6).
      number_of_iterations(30).
      number_of_relaxations(5).
      subdivision_ratio(0.5).
      face_proxy_map(fpxmap).
      proxies(std::back_inserter(proxies)).
      anchors(std::back_inserter(points)).
      triangles(std::back_inserter(triangles)));

  std::cout << "#proxies " << proxies.size() << std::endl;
  std::cout << "#vertices " << points.size() << std::endl;
  std::cout << "#triangles " << triangles.size() << std::endl;

  return EXIT_SUCCESS;
}
