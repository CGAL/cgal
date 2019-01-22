#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Variational_shape_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;

typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;

namespace PMP = CGAL::Polygon_mesh_processing;

bool check_strict_ordering(const std::vector<FT> &error)
{
  if (error.empty()) {
    std::cout << "Empty error sequence." << std::endl;
    return false;
  }
  FT pre = error.front();
  for (std::vector<FT>::const_iterator itr = error.begin(); itr != error.end(); ++itr)
    if (pre < *itr)
      return false;

  return true;
}

/**
 * This file tests the decrease of the relaxing error on a sphere shape.
 */
int main()
{
  Mesh mesh;
  std::ifstream input("./data/sphere.off");
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

  // algorithm instance
  L21_metric error_metric(mesh,
    get(boost::vertex_point, const_cast<Mesh &>(mesh)));
  L21_approx approx(mesh,
    get(boost::vertex_point, const_cast<Mesh &>(mesh)),
    error_metric);

  approx.initialize_seeds(CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::RANDOM)
    .max_number_of_proxies(100));
  std::vector<FT> error;
  for (std::size_t i = 0; i < 30; ++i) {
    approx.run();
    error.push_back(approx.compute_total_error());
  }

  if (check_strict_ordering(error)) {
    std::cout << "Pass the decrease test." << std::endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "Failed the decrease test." << std::endl;
    return EXIT_FAILURE;
  }
}
