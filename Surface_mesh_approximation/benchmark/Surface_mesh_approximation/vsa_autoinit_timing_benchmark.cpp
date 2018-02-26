#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/VSA_approximation.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::VSA_approximation<Polyhedron, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;
typedef L21_approx::Proxy_fitting L21_proxy_fitting;

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

  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }
  std::cerr << "#triangles " << mesh.size_of_facets() << std::endl;

  // algorithm instance
  L21_approx approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  // set metric error and fitting functors
  L21_metric error_metric(mesh);
  L21_proxy_fitting proxy_fitting(mesh);
  approx.set_metric(error_metric, proxy_fitting);

  int method = std::atoi(argv[2]);
  if (method < 0 || method > 2)
    return EXIT_FAILURE;
  const FT error_drop(std::atof(argv[3]));
  int nb_relaxations = std::atoi(argv[4]);
  std::cerr << "#method " << method << std::endl;
  std::cerr << "#error_drop " << error_drop << std::endl;
  std::cerr << "#nb_relaxations " << nb_relaxations << std::endl;

  Timer t;
  std::cerr << "start seeding" << std::endl;
  t.start();
  approx.seeding(
    static_cast<CGAL::Approximation_seeding_tag>(method), boost::none, error_drop, nb_relaxations);
  t.stop();
  std::cerr << "seeding time " << t.time() << " sec." << std::endl;
  std::cerr << "#proxies " << approx.proxies_size() << std::endl;

  return EXIT_SUCCESS;
}
