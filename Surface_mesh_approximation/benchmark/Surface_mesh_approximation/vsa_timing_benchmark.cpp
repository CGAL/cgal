#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/vsa_approximation.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::VSA::Mesh_approximation<Polyhedron, Vertex_point_map> L21_apporx;
typedef L21_apporx::Error_metric L21_metric;
typedef L21_apporx::Proxy_fitting L21_proxy_fitting;

typedef CGAL::Timer Timer;

/**
 * This file is a timing benchmark of each phase.
 * With different configuration:
   1. initialization
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

  // algorithm instance
  L21_apporx approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  // set metric error and fitting functors
  L21_metric error_metric(mesh);
  L21_proxy_fitting proxy_fitting(mesh);
  approx.set_metric(error_metric, proxy_fitting);

  int init = std::atoi(argv[2]);
  if (init < 0 || init > 2)
    return EXIT_FAILURE;
  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::cerr << "#init " << init << std::endl;
  std::cerr << "#num_proxies " << num_proxies << std::endl;
  std::cerr << "#num_iterations " << num_iterations << std::endl;

  Timer t0, t1;
  t1.start();

  std::cerr << "start initialization" << std::endl;
  t0.reset();
  t0.start();
  approx.init_by_number(
    static_cast<CGAL::VSA::Seeding>(init), num_proxies);
  t0.stop();
  std::cerr << "initialization time " << t0.time() << " sec." << std::endl;

  std::cerr << "start relaxation" << std::endl;
  t0.reset();
  t0.start();
  approx.run(num_iterations);
  t0.stop();
  std::cerr << "relaxation time " << t0.time() << " sec." << std::endl;

  Polyhedron mesh_out;
  t0.reset();
  t0.start();
  approx.extract_mesh(mesh_out);
  t0.stop();
  std::cerr << "meshing time " << t0.time() << " sec." << std::endl;

  std::cerr << "total time " << t1.time() << " sec." << std::endl;

  return EXIT_SUCCESS;
}
