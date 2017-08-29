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

bool test_manifold(const char *file_name, const FT drop = FT(1e-8))
{
  Polyhedron mesh;
  std::ifstream input(file_name);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cout << "Invalid off file." << std::endl;
    return false;
  }

  // algorithm instance
  L21VSA vsa_l21(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  L21Metric l21_metric(mesh);
  L21ProxyFitting l21_fitting(mesh);
  vsa_l21.set_metric(l21_metric, l21_fitting);

  // approximation, init from error, drop to the target error incrementally
  const std::size_t num_iterations = 20;
  const std::size_t inner_iterations = 5;
  vsa_l21.seeding_by_error(L21VSA::Incremental, drop, inner_iterations);
  for (std::size_t i = 0; i < num_iterations; ++i)
    vsa_l21.run_one_step();
  std::cout << "#proxies " << vsa_l21.get_proxies_size() << std::endl;

  // meshing
  Polyhedron mesh_out;
  return vsa_l21.meshing(mesh_out);
}

/**
 * This file tests the meshing of the algorithm.
 * For now, we can only expect manifold output on simple geometric objects.
 */
int main()
{
  std::cout << "Meshing manifold test." << std::endl;
  const char file0[] = "./data/cube_meshed.off";
  std::cout << "Testing file \"" << file0 << '\"' << std::endl;
  if (!test_manifold(file0)) {
    std::cout << "Failed." << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "Succeeded." << std::endl;

  const char file1[] = "./data/cube_meshed_open.off";
  std::cout << "Testing file \"" << file1 << '\"' << std::endl;
  if (!test_manifold(file1)) {
    std::cout << "Failed." << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "Succeeded." << std::endl;

  const char file2[] = "./data/sphere_iso.off";
  std::cout << "Testing file \"" << file2 << '\"' << std::endl;
  if (!test_manifold(file2, FT(1e-2))) {
    std::cout << "Failed." << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "Succeeded." << std::endl;

  return EXIT_SUCCESS;
}
