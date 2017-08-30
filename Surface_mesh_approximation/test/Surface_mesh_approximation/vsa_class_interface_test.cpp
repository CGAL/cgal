#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_handle Facet_handle;
typedef boost::associative_property_map<std::map<Facet_handle, std::size_t> > FacetProxyMap;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;

typedef CGAL::L2Metric<Polyhedron> L2Metric;
typedef CGAL::L2ProxyFitting<Polyhedron> L2ProxyFitting;
typedef CGAL::VSA_approximation<Polyhedron, VertexPointMap,
  L2Metric, L2ProxyFitting> L2VSA;
typedef L2VSA::Proxy PlaneProxy;

/**
 * This file tests the main class API and the L2 metric.
 * It should cover all the APIs.
 */
int main()
{
  Polyhedron mesh;
  std::ifstream input("./data/sphere_iso.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // facet area map
  std::map<Facet_handle, std::size_t> facet_index;
  for (Polyhedron::Facet_iterator fitr = mesh.facets_begin();
    fitr != mesh.facets_end(); ++fitr)
    facet_index.insert(std::pair<Facet_handle, std::size_t>(fitr, 0));
  FacetProxyMap proxy_pmap(facet_index);

  // create L2VSA L2 metric approximation algorithm instance
  std::cout << "setup algorithm instance" << std::endl;
  L2VSA l2_approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));
  L2Metric metric(mesh);
  L2ProxyFitting proxy_fitting(mesh);
  l2_approx.set_metric(metric, proxy_fitting);

  // random init and run
  std::cout << "random init and run" << std::endl;
  l2_approx.seeding_by_number(L2VSA::Random, 10);
  for (std::size_t i = 0; i < 10; ++i)
    l2_approx.run_one_step();
  if (l2_approx.get_proxies_size() != 10)
    return EXIT_FAILURE;

  // incremental add and run until convergence
  std::cout << "incremental add and run until convergence" << std::endl;
  l2_approx.add_proxies_furthest(3, 5);
  if (l2_approx.run_until_convergence(0.1))
    std::cout << "Converged." << std::endl;
  if (l2_approx.get_proxies_size() != 13)
    return EXIT_FAILURE;

  std::cout << "hierarchical add and run until convergence" << std::endl;
  l2_approx.add_proxies_error_diffusion(3);
  for (std::size_t i = 0; i < 10; ++i)
    l2_approx.run_one_step();
  if (l2_approx.get_proxies_size() != 16)
    return EXIT_FAILURE;

  // merge and teleport the proxies from local minimal
  std::cout << "teleport" << std::endl;
  l2_approx.teleport_proxies(3, false);
  for (std::size_t i = 0; i < 10; ++i)
    l2_approx.run_one_step();
  if (l2_approx.get_proxies_size() != 16)
    return EXIT_FAILURE;

  // split proxy 0 into 2 proxies
  // precondition: proxy 0 should have more than 2 facets
  std::cout << "spliting" << std::endl;
  l2_approx.split(0);
  for (std::size_t i = 0; i < 10; ++i)
    l2_approx.run_one_step();
  if (l2_approx.get_proxies_size() != 17)
    return EXIT_FAILURE;

  // extract the approximation polyhedron
  std::cout << "meshing" << std::endl;
  Polyhedron out_mesh;
  if (l2_approx.meshing(out_mesh, FT(0.5), true))
    std::cout << "manifold." << std::endl;
  else
    std::cout << "non-manifold" << std::endl;

  // get outputs
  std::cout << "get outputs" << std::endl;
  l2_approx.get_proxy_map(proxy_pmap);

  std::vector<PlaneProxy> proxies;
  l2_approx.get_proxies(std::back_inserter(proxies));

  std::vector<Point> anchor_pos;
  l2_approx.get_anchor_points(std::back_inserter(anchor_pos));

  std::vector<Polyhedron::Vertex_handle> anchor_vtx;
  l2_approx.get_anchor_vertices(std::back_inserter(anchor_vtx));

  std::vector<int> tris;
  l2_approx.get_indexed_triangles(std::back_inserter(tris));

  std::vector<std::vector<std::size_t> > boundary;
  boundary = l2_approx.get_indexed_boundary_polygons();

  const FT drop(0.001);
  const std::size_t iterations = 5;
  std::cout << "rebuild and hierarchical init" << std::endl;
  l2_approx.rebuild();
  if (l2_approx.get_proxies_size() != 0)
    return EXIT_FAILURE;
  l2_approx.seeding_by_error(L2VSA::Hierarchical, drop, iterations);
  for (std::size_t i = 0; i < 10; ++i)
    l2_approx.run_one_step();
  std::cout << "#proxies " << l2_approx.get_proxies_size() << std::endl;

  std::cout << "rebuild and incremental init" << std::endl;
  l2_approx.rebuild();
  if (l2_approx.get_proxies_size() != 0)
    return EXIT_FAILURE;
  l2_approx.seeding_by_error(L2VSA::Incremental, drop, iterations);
  for (std::size_t i = 0; i < 10; ++i)
    l2_approx.run_one_step();
  std::cout << "#proxies " << l2_approx.get_proxies_size() << std::endl;

  return EXIT_SUCCESS;
}
