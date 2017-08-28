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
    std::cerr << "Invalid off file." << std::endl;
    return false;
  }

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
  vsa_l21.seeding_by_error(L21VSA::IncrementalInit, drop, inner_iterations);
  for (std::size_t i = 0; i < num_iterations; ++i)
    vsa_l21.run_one_step();
  if (vsa_l21.get_proxies_size() != target_num_proxies) {
    std::cout << "#targeted - #result "
      << target_num_proxies << ' '
      << vsa_l21.get_proxies_size() << std::endl;
    std::cout << "incremental reaching failed" << std::endl;
    return false;
  }

  // meshing, should be manifold
  Polyhedron mesh_out;
  if (!vsa_l21.meshing(mesh_out)) {
    std::cout << "incremental reaching meshing non-manifold" << std::endl;
    return false;
  }

  return true;
}

/**
 * This file tests the correctness of the algorithm.
 * Basically we input a cube mesh and see if it outputs a cube.
 */
int main()
{
  const char file0[] = "./data/cube_meshed.off";
  std::cout << "Testing file \"" << file0 << '\"';
  if (!test_shape(file0, 6)) {
    std::cout << "Failed." << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "Succeeded." << std::endl;

  const char file1[] = "./data/cube_meshed_open.off";
  std::cout << "Testing file \"" << file1 << '\"';
  if (!test_shape(file1, 5)) {
    std::cout << "Failed." << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "Succeeded." << std::endl;

  return EXIT_SUCCESS;
}
