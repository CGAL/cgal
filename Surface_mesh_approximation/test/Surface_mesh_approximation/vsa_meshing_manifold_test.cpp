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
  L21_approx approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  L21_metric error_metric(mesh);
  L21_proxy_fitting proxy_fitting(mesh);
  approx.set_metric(error_metric, proxy_fitting);

  // approximation, seeding from error, drop to the target error incrementally
  const std::size_t num_iterations = 20;
  const std::size_t inner_iterations = 5;
  approx.seeding(CGAL::VSA::Incremental, boost::none, drop, inner_iterations);
  approx.run(num_iterations);
  std::cout << "#proxies " << approx.proxies_size() << std::endl;

  // meshing
  Polyhedron mesh_out;
  if (approx.extract_mesh(mesh_out)) {
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
