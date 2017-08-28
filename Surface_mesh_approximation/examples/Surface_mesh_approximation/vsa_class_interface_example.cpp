#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;

typedef CGAL::VSA_approximation<Polyhedron, VertexPointMap> L21VSA;
typedef L21VSA::ErrorMetric L21Metric;
typedef L21VSA::ProxyFitting L21ProxyFitting;

int main()
{
  // read Polyhedron
  Polyhedron input;
  std::ifstream file("data/bear.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // create VSA L21 metric approximation algorithm instance
  L21VSA l21_approx(input,
    get(boost::vertex_point, const_cast<Polyhedron &>(input)));
  // set error and fitting functors
  L21Metric metric(input);
  L21ProxyFitting proxy_fitting(input);
  l21_approx.set_metric(metric, proxy_fitting);

  // initialize proxies randomly on the mesh
  l21_approx.seeding_by_number(L21VSA::Random, 100);
  
  // run the iteration to minimize the error
  for (std::size_t i = 0; i < 30; ++i)
    l21_approx.run_one_step();

  // add proxies to the one with the maximum fitting error
  l21_approx.add_proxies_furthest(3);
  for (std::size_t i = 0; i < 10; ++i)
    l21_approx.run_one_step();

  // merge and teleport the proxies from local minimal
  l21_approx.teleport_proxies(2);
  for (std::size_t i = 0; i < 10; ++i)
    l21_approx.run_one_step();

  // extract the approximation polyhedron
  Polyhedron output;
  l21_approx.meshing(output);

  return EXIT_SUCCESS;
}
