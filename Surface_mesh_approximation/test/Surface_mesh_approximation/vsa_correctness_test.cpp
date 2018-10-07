#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Variational_shape_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;
typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;

namespace PMP = CGAL::Polygon_mesh_processing;

bool load_mesh(const char *file_name, Mesh &mesh)
{
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
  std::cout << "Load mesh " << file_name << " done." << std::endl;

  return true;
}

bool test_shape(const Mesh &mesh, const std::size_t target_num_proxies)
{
  // algorithm instance
  L21_metric error_metric(mesh,
    get(boost::vertex_point, const_cast<Mesh &>(mesh)));
  L21_approx approx(mesh,
    get(boost::vertex_point, const_cast<Mesh &>(mesh)),
    error_metric);

  // approximation, seeding from error, drop to the target error incrementally
  // should reach targeted number of proxies gradually
  const Kernel::FT drop(1e-2);
  const std::size_t num_iterations = 20;
  const std::size_t inner_iterations = 10;
  approx.initialize_seeds(
    CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::INCREMENTAL)
    .min_error_drop(drop)
    .number_of_relaxations(inner_iterations));
  approx.run(num_iterations);

  // eliminate redundant area (local minima) by merging
  boost::optional<std::pair<std::size_t, std::size_t> > best_pair = boost::none;
  while ((best_pair = approx.find_best_merge(true)) != boost::none) {
    approx.merge(best_pair->first, best_pair->second);
    approx.run(num_iterations);
  }

  if (approx.number_of_proxies() != target_num_proxies) {
    std::cout << "#targeted - #result "
      << target_num_proxies << ' '
      << approx.number_of_proxies() << std::endl;

    std::cout << "Failed." << std::endl;
    return false;
  }

  std::cout << "Succeeded." << std::endl;
  return true;
}

/**
 * This file tests the correctness of the algorithm.
 * The correctness is verified by seeding by error
 * and check if the number of desired proxies are generated.
 * Basically we input a cube mesh and see if it outputs 6 proxies.
 */
int main()
{
  const char file_cube[] = "./data/cube.off";
  std::cout << "Testing close mesh " << file_cube << std::endl;
  Mesh mesh_cube;
  if (!load_mesh(file_cube, mesh_cube) || !test_shape(mesh_cube, 6))
    return EXIT_FAILURE;

  const char file_cube2[] = "./data/cube-ouvert.off";
  std::cout << "Testing open mesh " << file_cube2 << std::endl;
  Mesh mesh_cube2;
  if (!load_mesh(file_cube2, mesh_cube2) || !test_shape(mesh_cube2, 5))
    return EXIT_FAILURE;

  std::cout << "Tesh mesh with disconnected components" << std::endl;
  Mesh mesh_merged = mesh_cube;
  mesh_cube2.collect_garbage();
  // the second parameter of operator+= should not have garbage, or merge will crash
  mesh_merged += mesh_cube2;
  std::cout << "Mege done \n#F "
    << std::distance(faces(mesh_merged).first, faces(mesh_merged).second)
    << "\n#V " << std::distance(vertices(mesh_merged).first, vertices(mesh_merged).second)
    <<  std::endl;
  if (!test_shape(mesh_merged, 11))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
