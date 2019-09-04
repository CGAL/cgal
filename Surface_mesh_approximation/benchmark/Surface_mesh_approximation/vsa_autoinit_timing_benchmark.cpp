#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;

typedef CGAL::Timer Timer;

/**
 * This file is a timing benchmark of the automatic seeding.
 * With different configuration:
   1. seeding method
   2. error drop tolerance
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
  const FT error_drop(std::atof(argv[3]));
  int number_of_relaxations = std::atoi(argv[4]);
  std::cout << "#method " << method << std::endl;
  std::cout << "#error_drop " << error_drop << std::endl;
  std::cout << "#number_of_relaxations " << number_of_relaxations << std::endl;

  Timer t;
  std::cout << "start seeding" << std::endl;
  t.start();
  approx.initialize_seeds(
    CGAL::parameters::seeding_method(static_cast<CGAL::Surface_mesh_approximation::Seeding_method>(method))
    .min_error_drop(error_drop)
    .number_of_relaxations(number_of_relaxations));
  t.stop();
  std::cout << "seeding time " << t.time() << " sec." << std::endl;
  std::cout << "#proxies " << approx.number_of_proxies() << std::endl;

  return EXIT_SUCCESS;
}
