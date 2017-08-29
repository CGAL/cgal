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

  std::cout << "Testing \"" << file_name << '\"' << std::endl;
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
  if (vsa_l21.meshing(mesh_out)) {
    std::cout << "Succeeded." << std::endl;
    return true;
  }

  std::cout << "Failed." << std::endl;
  return false;
}

/**
 * This file tests the meshing of the algorithm.
 * For now, we can only expect manifold output on simple geometric objects.
 */
int main()
{
  std::cout << "Meshing manifold test." << std::endl;
  if (!test_manifold("./data/cube_meshed.off"))
    return EXIT_FAILURE;

  if (!test_manifold("./data/cube_meshed_open.off"))
    return EXIT_FAILURE;

  if (!test_manifold("./data/sphere_iso.off", FT(1e-2)))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
