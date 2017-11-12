#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>
#include <CGAL/vsa_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Point_3 Point;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;
typedef boost::associative_property_map<std::map<Facet_handle, Point> > FacetCenterMap;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;

// user defined point-wise compact metric
struct CompactMetric {
  typedef Point Proxy;

  CompactMetric(const FacetCenterMap &_center_pmap)
    : center_pmap(_center_pmap) {}

  FT operator()(const Facet_handle &f, const Proxy &px) const {
    return FT(std::sqrt(CGAL::to_double(
      CGAL::squared_distance(center_pmap[f], px))));
  }

  const FacetCenterMap center_pmap;
};

struct PointProxyFitting {
  typedef Point Proxy;

  PointProxyFitting(const FacetCenterMap &_center_pmap, const FacetAreaMap &_area_pmap)
    : center_pmap(_center_pmap), area_pmap(_area_pmap) {}

  template<typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    // fitting center
    Vector center = CGAL::NULL_VECTOR;
    FT area(0);
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      center = center + (center_pmap[*fitr] - CGAL::ORIGIN) * area_pmap[*fitr];
      area += area_pmap[*fitr];
    }
    center = center / area;
    return CGAL::ORIGIN + center;
  }

  const FacetCenterMap center_pmap;
  const FacetAreaMap area_pmap;
};
typedef CGAL::VSA_approximation<Polyhedron, VertexPointMap,
  CompactMetric, PointProxyFitting> CompactVSA;

/**
 * This file tests the user defined metric.
 */
int main()
{
  Polyhedron mesh;
  std::ifstream input("./data/cube_meshed_open.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // construct facet normal & area map
  std::map<Facet_handle, FT> facet_areas;
  std::map<Facet_handle, Point> facet_centers;
  for(Facet_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    const Halfedge_handle he = fitr->halfedge();
    const Point &p0 = he->opposite()->vertex()->point();
    const Point &p1 = he->vertex()->point();
    const Point &p2 = he->next()->vertex()->point();
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<Facet_handle, FT>(fitr, area));
    facet_centers.insert(std::pair<Facet_handle, Point>(fitr, CGAL::centroid(p0, p1, p2)));
  }
  FacetAreaMap area_pmap(facet_areas);
  FacetCenterMap center_pmap(facet_centers);

  // create compact metric approximation algorithm instance
  std::cout << "create compact vas instance" << std::endl;
  CompactVSA compact_approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));

  CompactMetric metric(center_pmap);
  PointProxyFitting proxy_fitting(center_pmap, area_pmap);
  compact_approx.set_metric(metric, proxy_fitting);

  std::cout << "random init and run" << std::endl;
  compact_approx.init_by_number(CGAL::VSA_seeding::Random, 20);
  for (std::size_t i = 0; i < 20; ++i)
    compact_approx.run_one_step();
  if (compact_approx.get_proxies_size() != 20)
    return EXIT_FAILURE;

  // extract the approximation polyhedron
  std::cout << "meshing" << std::endl;
  Polyhedron out_mesh;
  if (compact_approx.meshing(out_mesh))
    std::cout << "manifold." << std::endl;
  else
    std::cout << "non-manifold" << std::endl;

  return EXIT_SUCCESS;
}
