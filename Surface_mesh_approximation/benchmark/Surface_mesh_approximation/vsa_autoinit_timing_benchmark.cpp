#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/vsa_approximation.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::VSA::Mesh_approximation<Polyhedron, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;
typedef L21_approx::Proxy_fitting L21_proxy_fitting;

typedef CGAL::Timer Timer;

/**
 * This file is a timing benchmark of the automatic initialization.
 * With different configuration:
   1. initialization
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

  int init = std::atoi(argv[2]);
  if (init < 0 || init > 2)
    return EXIT_FAILURE;
  const FT tol(std::atof(argv[3]));
  int iterations = std::atoi(argv[4]);
  std::cerr << "#init " << init << std::endl;
  std::cerr << "#tolerance " << tol << std::endl;
  std::cerr << "#iterations " << iterations << std::endl;

  Timer t;
  std::cerr << "start initialization" << std::endl;
  t.start();
  approx.init_by_error(
    static_cast<CGAL::VSA::Seeding>(init), tol, iterations);
  t.stop();
  std::cerr << "initialization time " << t.time() << " sec." << std::endl;
  std::cerr << "#proxies " << approx.get_proxies_size() << std::endl;

  return EXIT_SUCCESS;
}
