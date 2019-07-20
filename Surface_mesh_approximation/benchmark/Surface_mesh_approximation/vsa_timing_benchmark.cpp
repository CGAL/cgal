#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;

typedef CGAL::Timer Timer;

/**
 * This file is a timing benchmark of each phase.
 * With different configuration:
   1. seeding method
   2. number of proxies
   3. number of iterations
 * TODO: error
 */
int main(int argc, char *argv[])
{
  if (argc < 5)
    return EXIT_FAILURE;

  Mesh mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cout << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "#faces "
    << std::distance(faces(mesh).first, faces(mesh).second) << std::endl;

  Vertex_point_map vpmap = get(boost::vertex_point, const_cast<Mesh &>(mesh));

  // error metric and fitting functors
  L21_metric error_metric(mesh, vpmap);

  // algorithm instance
  L21_approx approx(mesh, vpmap, error_metric);

  int method = std::atoi(argv[2]);
  if (method < 0 || method > 2)
    return EXIT_FAILURE;
  const std::size_t nb_proxies = std::atoi(argv[3]);
  const std::size_t nb_iterations = std::atoi(argv[4]);
  std::cout << "#method " << method << std::endl;
  std::cout << "#nb_proxies " << nb_proxies << std::endl;
  std::cout << "#nb_iterations " << nb_iterations << std::endl;

  Timer t0, t1;
  t1.start();

  std::cout << "start seeding" << std::endl;
  t0.reset();
  t0.start();
  approx.initialize_seeds(
    CGAL::parameters::seeding_method(static_cast<CGAL::Surface_mesh_approximation::Seeding_method>(method))
    .max_number_of_proxies(nb_proxies));
  t0.stop();
  std::cout << "seeding time " << t0.time() << " sec." << std::endl;

  std::cout << "start iterations" << std::endl;
  t0.reset();
  t0.start();
  approx.run(nb_iterations);
  t0.stop();
  std::cout << "iterations time " << t0.time() << " sec." << std::endl;

  t0.reset();
  t0.start();
  approx.extract_mesh(CGAL::parameters::all_default());
  t0.stop();
  std::cout << "meshing time " << t0.time() << " sec." << std::endl;

  std::cout << "total time " << t1.time() << " sec." << std::endl;

  return EXIT_SUCCESS;
}
