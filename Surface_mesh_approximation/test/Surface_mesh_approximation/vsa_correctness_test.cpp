#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/vsa_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::VSA::Mesh_approximation<Polyhedron, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;
typedef L21_approx::Proxy_fitting L21_proxy_fitting;

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
  L21_approx approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  L21_metric error_metric(mesh);
  L21_proxy_fitting proxy_fitting(mesh);
  approx.set_metric(error_metric, proxy_fitting);

  // approximation, init from error, drop to the target error incrementally
  // should reach targeted number of proxies gradually
  const FT drop(1e-8);
  const std::size_t num_iterations = 20;
  const std::size_t inner_iterations = 10;
  approx.init_by_error(CGAL::VSA::Incremental, drop, inner_iterations);
  approx.run(num_iterations);

  // eliminate redundant area (local minima) by merging
  std::size_t px0 = 0, px1 = 0;
  while (approx.find_best_merge(px0, px1, true)) {
    approx.merge(px0, px1);
    approx.run(num_iterations);
  }

  if (approx.get_proxies_size() != target_num_proxies) {
    std::cout << "#targeted - #result "
      << target_num_proxies << ' '
      << approx.get_proxies_size() << std::endl;

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
