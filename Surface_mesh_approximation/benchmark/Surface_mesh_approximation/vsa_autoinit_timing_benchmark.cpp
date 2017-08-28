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
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;

typedef CGAL::VSA_approximation<Polyhedron, VertexPointMap> L21VSA;
typedef L21VSA::ErrorMetric L21Metric;
typedef L21VSA::ProxyFitting L21ProxyFitting;

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
    return 1;

  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return 1;
  }
  std::cerr << "#triangles " << mesh.size_of_facets() << std::endl;

  // algorithm instance
  L21VSA l21_vsa(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  // set metric error and fitting functors
  L21Metric l21_metric(mesh);
  L21ProxyFitting l21_fitting(mesh);
  l21_vsa.set_metric(l21_metric, l21_fitting);

  int init = std::atoi(argv[2]);
  if (init < 0 || init > 2)
    return 1;
  const FT tol(std::atof(argv[3]));
  int iterations = std::atoi(argv[4]);
  std::cerr << "#init " << init << std::endl;
  std::cerr << "#tolerance " << tol << std::endl;
  std::cerr << "#iterations " << iterations << std::endl;

  Timer t;
  std::cerr << "start initialization" << std::endl;
  t.start();
  l21_vsa.seeding_error(
    static_cast<L21VSA::Initialization>(init),
    tol,
    iterations);
  t.stop();
  std::cerr << "initialization time " << t.time() << " sec." << std::endl;
  std::cerr << "#proxies " << l21_vsa.get_proxies_size() << std::endl;

  return 0;
}
