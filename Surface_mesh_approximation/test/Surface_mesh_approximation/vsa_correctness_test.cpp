#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;

typedef CGAL::VSA_approximation<Polyhedron, VertexPointMap> L21VSA;
typedef L21VSA::ErrorMetric L21Metric;
typedef L21VSA::ProxyFitting L21ProxyFitting;

bool test_shape(const char *file_name, const std::size_t target_num_proxies)
{
  Polyhedron mesh;
  std::ifstream input(file_name);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cout << "Invalid off file." << std::endl;
    return false;
  }

  std::cout << "Testing \"" << file_name << '\"' << std::endl;
  // algorithm instance
  L21VSA vsa_l21(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  L21Metric l21_metric(mesh);
  L21ProxyFitting l21_fitting(mesh);
  vsa_l21.set_metric(l21_metric, l21_fitting);

  // approximation, init from error, drop to the target error incrementally
  // should reach targeted number of proxies gradually
  const FT drop(1e-8);
  const std::size_t num_iterations = 20;
  const std::size_t inner_iterations = 10;
  vsa_l21.seeding_by_error(L21VSA::Incremental, drop, inner_iterations);
  for (std::size_t i = 0; i < num_iterations; ++i)
    vsa_l21.run_one_step();

  // eliminate redundant area (local minima) by merging
  std::size_t px0 = 0, px1 = 0;
  while (vsa_l21.find_best_merge(px0, px1, true)) {
    vsa_l21.merge(px0, px1);
    for (std::size_t i = 0; i < num_iterations; ++i)
      vsa_l21.run_one_step();
  }

  if (vsa_l21.get_proxies_size() != target_num_proxies) {
    std::cout << "#targeted - #result "
      << target_num_proxies << ' '
      << vsa_l21.get_proxies_size() << std::endl;

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
  std::cout << "Correctness test." << std::endl;
  if (!test_shape("./data/cube_meshed.off", 6))
    return EXIT_FAILURE;

  if (!test_shape("./data/cube_meshed_open.off", 5))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
