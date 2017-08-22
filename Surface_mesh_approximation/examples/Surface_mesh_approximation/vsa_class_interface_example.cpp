#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/vsa_mesh_approximation_metrics.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;

typedef CGAL::PlaneProxy<Polyhedron_3> PlaneProxy;
typedef CGAL::L2Metric<Polyhedron_3> L2Metric;
typedef CGAL::L2ProxyFitting<Polyhedron_3> L2ProxyFitting;
typedef CGAL::VSA_approximation<Polyhedron_3, PlaneProxy, L2Metric, L2ProxyFitting> VSA;

int main()
{
  // read Polyhedron_3
  Polyhedron_3 mesh;
  std::ifstream input("data/bear.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  L2Metric metric(mesh);
  L2ProxyFitting proxy_fitting(mesh);

  // create VSA L2 metric approximation algorithm instance
  VSA l2_approx;
  l2_approx.set_mesh(mesh);
  l2_approx.set_error_metric(metric);
  l2_approx.set_proxy_fitting(proxy_fitting);

  // initialize proxies randomly on the mesh
  l2_approx.init_proxies(100, VSA::RandomInit);
  
  // run the iteration to minimize the error
  for (std::size_t i = 0; i < 30; ++i)
    l2_approx.run_one_step();

  // add proxies to the one with the maximum fitting error
  l2_approx.add_proxies(VSA::IncrementalInit, 3);
  for (std::size_t i = 0; i < 10; ++i)
    l2_approx.run_one_step();

  // merge and teleport the proxies from local minimal
  l2_approx.teleport_proxies(2);
  for (std::size_t i = 0; i < 10; ++i)
    l2_approx.run_one_step();

  // extract the approximation polyhedron
  Polyhedron_3 out_mesh;
  l2_approx.meshing(out_mesh);

  return EXIT_SUCCESS;
}
