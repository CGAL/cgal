#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Halfedge_handle Halfedge_handle;

typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;
typedef boost::associative_property_map<std::map<Facet_handle, std::size_t> > FacetProxyMap;

typedef CGAL::VSA_approximation<Polyhedron, VertexPointMap> L21VSA;
typedef L21VSA::ErrorMetric L21Metric;
typedef L21VSA::ProxyFitting L21ProxyFitting;
typedef L21VSA::Proxy L21Proxy;

#define CGAL_VSA_TEST_TOLERANCE 1e-8

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
 * from planar part to the spherical one and lower the error.
 */
int main()
{
  Polyhedron mesh;
  std::ifstream input("./data/plane-sphere-high.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Teleportation test." << std::endl;

  // algorithm instance
  L21Metric l21_metric(mesh);
  L21ProxyFitting l21_fitting(mesh);
  L21VSA l21_vsa(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));
  l21_vsa.set_metric(l21_metric, l21_fitting);

  std::cout << "Seeding by number." << std::endl;
  l21_vsa.init_by_number(CGAL::VSA_seeding::Random, 50);
  if (l21_vsa.get_proxies_size() != 50)
    return EXIT_FAILURE;
  for (std::size_t i = 0; i < 10; ++i) {
    l21_vsa.partition();
    l21_vsa.fit();
  }

  // teleport until merge test failed
  std::vector<FT> error;
  std::size_t count = 0;
  while(l21_vsa.teleport_proxies(1) == 1) {
    FT sum_err(0);
    for (std::size_t i = 0; i < 10; ++i)
      sum_err += l21_vsa.run_one_step();
    error.push_back(sum_err / FT(10));
    ++count;
  }
  std::cout << "#teleportation " << count << std::endl;

  if (!check_strict_ordering(error)) {
    std::cout << "Failed: teleportation error decrease inconsistent." << std::endl;
    return EXIT_FAILURE;
  }

  std::map<Facet_handle, std::size_t> internal_fidxmap;
  for (Facet_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr)
    internal_fidxmap[fitr] = 0;
  FacetProxyMap fproxymap(internal_fidxmap);
  l21_vsa.get_proxy_map(fproxymap);
  std::vector<L21Proxy> proxies;
  l21_vsa.get_proxies(std::back_inserter(proxies));

  CGAL::Bbox_3 bbox = CGAL::bbox_3(mesh.points_begin(), mesh.points_end());
  const FT ymin = bbox.ymin(), ymax = bbox.ymax(), yrange = ymax - ymin;
  std::cout << "Range along y axis: [" << ymin << ", " << ymax << "]" << std::endl;

  // test if all facets on the planar part are in the same proxy
  std::size_t planar_pxidx = static_cast<std::size_t>(-1);
  std::size_t num_planar_facets = 0;
  bool first = true;
  for (Facet_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    Halfedge_handle he = fitr->halfedge();
    const Point &p0 = he->opposite()->vertex()->point();
    const Point &p1 = he->vertex()->point();
    const Point &p2 = he->next()->vertex()->point();
    const Point fcenter = CGAL::centroid(p0, p1, p2);
    Vector fnormal = CGAL::normal(p0, p1, p2);
    fnormal = fnormal / FT(std::sqrt(CGAL::to_double(fnormal.squared_length())));

    // check the facet center and normal to see if it is on the planar part of the geometry
    double dis_var = std::abs(CGAL::to_double((fcenter.y() - ymin) / yrange));
    double dir_var = std::abs(CGAL::to_double(fnormal.y()) - 1.0);
    if (dis_var < CGAL_VSA_TEST_TOLERANCE && dir_var < CGAL_VSA_TEST_TOLERANCE) {
      ++num_planar_facets;
      const std::size_t pxidx = fproxymap[fitr];
      if (first) {
        first = false;
        planar_pxidx = pxidx;
      }

      // only one proxy on the planar part of the geometry
      if (pxidx != planar_pxidx) {
        std::cout << "Failed: more than one proxy on the planar part of the model." << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  std::cout << "#facets on planar part " << num_planar_facets << std::endl;

  // proxy of the planar part should facing straight towards the y positive.
  double px_dir_var = std::abs(CGAL::to_double(proxies[planar_pxidx].normal.y()) - 1.0);
  if (px_dir_var > CGAL_VSA_TEST_TOLERANCE) {
    std::cout << "Failed: the proxy of planar part is incorrect." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Succeeded." << std::endl;
  return EXIT_SUCCESS;
}
