#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/vsa_mesh_approximation_metrics.h>
#include <CGAL/VSA_approximation.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;

typedef CGAL::PlaneProxy<Polyhedron_3> PlaneProxy;
typedef CGAL::L21Metric<Polyhedron_3> L21Metric;
typedef CGAL::L21ProxyFitting<Polyhedron_3> L21ProxyFitting;
typedef CGAL::VSA_approximation<Polyhedron_3, PlaneProxy, L21Metric, L21ProxyFitting> VSAL21;

typedef CGAL::Timer Timer;

/**
 * This file is a timing benchmark of the automatic initialization.
 * With different configuration:
   1. initialization
   2. error drop tolerance
 */
int main(int argc, char *argv[])
{
  if (argc < 4)
    return 1;

  Polyhedron_3 mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return 1;
  }
  std::cerr << "#triangles " << mesh.size_of_facets() << std::endl;

  L21Metric l21_metric(mesh);
  L21ProxyFitting l21_fitting(mesh);

  // algorithm instance
  VSAL21 vsa_l21(l21_metric, l21_fitting);
  vsa_l21.set_mesh(mesh);

  int init = std::atoi(argv[2]);
  if (init < 0 || init > 2)
    return 1;
  const FT tol(std::atof(argv[3]));
  std::cerr << "#init " << init << std::endl;
  std::cerr << "#tolerance " << tol << std::endl;

  Timer t;
  std::cerr << "start initialization" << std::endl;
  t.start();
  vsa_l21.init_proxies_error(tol, static_cast<VSAL21::Initialization>(init));
  t.stop();
  std::cerr << "initialization time " << t.time() << " sec." << std::endl;
  std::cerr << "#proxies " << vsa_l21.get_proxies_size() << std::endl;

  return 0;
}
