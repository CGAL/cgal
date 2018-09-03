#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::Variational_shape_approximation<Polyhedron, Vertex_point_map> L21_approx;
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

  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }
  std::cerr << "#triangles " << mesh.size_of_facets() << std::endl;

  // error metric and fitting functors
  L21_metric error_metric(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  // algorithm instance
  L21_approx approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)),
    error_metric);

  int method = std::atoi(argv[2]);
  if (method < 0 || method > 2)
    return EXIT_FAILURE;
  const std::size_t nb_proxies = std::atoi(argv[3]);
  const std::size_t nb_iterations = std::atoi(argv[4]);
  std::cerr << "#method " << method << std::endl;
  std::cerr << "#nb_proxies " << nb_proxies << std::endl;
  std::cerr << "#nb_iterations " << nb_iterations << std::endl;

  Timer t0, t1;
  t1.start();

  std::cerr << "start seeding" << std::endl;
  t0.reset();
  t0.start();
  approx.initialize_seeds(
    static_cast<CGAL::VSA::Seeding_method>(method), nb_proxies);
  t0.stop();
  std::cerr << "seeding time " << t0.time() << " sec." << std::endl;

  std::cerr << "start iterations" << std::endl;
  t0.reset();
  t0.start();
  approx.run(nb_iterations);
  t0.stop();
  std::cerr << "iterations time " << t0.time() << " sec." << std::endl;

  t0.reset();
  t0.start();
  approx.extract_mesh(CGAL::VSA::parameters::all_default());
  t0.stop();
  std::cerr << "meshing time " << t0.time() << " sec." << std::endl;

  std::cerr << "total time " << t1.time() << " sec." << std::endl;

  return EXIT_SUCCESS;
}
