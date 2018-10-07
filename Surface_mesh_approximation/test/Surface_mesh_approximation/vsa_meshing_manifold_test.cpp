#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Variational_shape_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;

namespace PMP = CGAL::Polygon_mesh_processing;

bool test_manifold(const char *file_name, const FT drop = FT(1e-2))
{
  Mesh mesh;
  std::ifstream input(file_name);
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cout << "Invalid input file." << std::endl;
    return false;
  }

  const std::size_t nb_removed = PMP::remove_isolated_vertices(mesh);
  if (nb_removed > 0)
    std::cout << nb_removed << " isolated vertices are removed." << std::endl;

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

  std::cout << "Testing \"" << file_name << '\"' << std::endl;
  // algorithm instance
  L21_metric error_metric(mesh,
    get(boost::vertex_point, const_cast<Mesh &>(mesh)));
  L21_approx approx(mesh,
    get(boost::vertex_point, const_cast<Mesh &>(mesh)),
    error_metric);

  // approximation, seeding from error, drop to the target error incrementally
  const std::size_t num_iterations = 20;
  const std::size_t inner_iterations = 5;
  approx.initialize_seeds(
    CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::INCREMENTAL)
    .min_error_drop(drop)
    .number_of_relaxations(inner_iterations));
  approx.run(num_iterations);
  std::cout << "#proxies " << approx.number_of_proxies() << std::endl;

  // meshing
  if (approx.extract_mesh(CGAL::parameters::subdivision_ratio(5.0))) {
    std::cout << "Succeeded." << std::endl;
    return true;
  }

  std::cout << "Failed." << std::endl;
  return false;
}

/**
 * This file tests the meshing of the algorithm.
 * For now, we can only expect manifold output on simple geometric objects.
 */
int main()
{
  std::cout << "Meshing manifold test." << std::endl;
  if (!test_manifold("./data/cube.off"))
    return EXIT_FAILURE;

  if (!test_manifold("./data/cube-ouvert.off"))
    return EXIT_FAILURE;

  if (!test_manifold("./data/sphere.off"))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
