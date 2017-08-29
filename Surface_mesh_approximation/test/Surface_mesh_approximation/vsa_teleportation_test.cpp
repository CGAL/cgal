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
typedef L21VSA::Proxy Proxy;

bool check_strict_ordering(const std::vector<FT> &error)
{
  if (error.empty()) {
    std::cout << "Empty error sequence." << std::endl;
    return false;
  }
  FT pre = error.front();
  for (std::vector<FT>::const_iterator itr = error.begin(); itr != error.end(); ++itr)
    if (pre < *itr)
      return false;

  return true;
}

/**
 * This file tests the teleportation functionality on a plane-sphere shape.
 * We initialize random first, then verify that teleporting all (most) planes
 * from planar parts to the spherical one and lower the error.
 * TODO: check if all planar parts to the spherical one
 */
int main()
{
  Polyhedron mesh;
  std::ifstream input("./data/plane-sphere-high.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // algorithm instance
  L21Metric l21_metric(mesh);
  L21ProxyFitting l21_fitting(mesh);
  L21VSA l21_vsa(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));
  l21_vsa.set_metric(l21_metric, l21_fitting);

  std::cout << "Seeding by number." << std::endl;
  l21_vsa.seeding_by_number(L21VSA::Random, 50);
  if (l21_vsa.get_proxies_size() != 50)
    return EXIT_FAILURE;
  for (std::size_t i = 0; i < 10; ++i) {
    l21_vsa.partition();
    l21_vsa.fit();
  }

  // teleport until merge test failed
  std::vector<FT> error;
  std::size_t count = 0;
  while(l21_vsa.teleport_proxies(1, true) == 1) {
    FT sum_err(0);
    for (std::size_t i = 0; i < 10; ++i)
      sum_err += l21_vsa.run_one_step();
    error.push_back(sum_err / FT(10));
    ++count;
  }
  std::cout << "#teleportation " << count << std::endl;

  if (check_strict_ordering(error)) {
    std::cout << "Pass the teleportation decrease test." << std::endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "Failed the teleportation decrease test." << std::endl;
    return EXIT_FAILURE;
  }
}
