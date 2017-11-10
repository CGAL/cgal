#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Point_3 Point;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Facet_iterator Facet_iterator;

typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;
typedef boost::associative_property_map<std::map<Facet_handle, Point> > FacetCenterMap;

// use a point as proxy
typedef Point PointProxy;

// user defined point-wise compact metric
struct CompactMetric {
  typedef PointProxy Proxy;

  // keeping a precomputed property map to speed up the computation
  CompactMetric(const FacetCenterMap &_center_pmap)
    : center_pmap(_center_pmap) {}

  // calculate the error of a facet to a proxy
  // here is just the Euclidean distance
  FT operator()(const Facet_handle &f, const Proxy &px) const {
    return FT(std::sqrt(CGAL::to_double(
      CGAL::squared_distance(center_pmap[f], px))));
  }

  const FacetCenterMap center_pmap;
};

// proxy fitting functor
struct PointProxyFitting {
  typedef PointProxy Proxy;

  // keeping precomputed property maps to speed up the computation
  PointProxyFitting(const FacetCenterMap &_center_pmap, const FacetAreaMap &_area_pmap)
    : center_pmap(_center_pmap), area_pmap(_area_pmap) {}

  // a template functor fit a new proxy from a range of facets
  template<typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    // fitting center
    Vector center = CGAL::NULL_VECTOR;
    FT area = FT(0.0);
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

int main()
{
  // create and read Polyhedron
  Polyhedron input;
  std::ifstream file("data/bear.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // construct precomputed facet normal and area map
  std::map<Facet_handle, FT> facet_areas;
  std::map<Facet_handle, Point> facet_centers;
  for(Facet_iterator fitr = input.facets_begin(); fitr != input.facets_end(); ++fitr) {
    const Halfedge_handle he = fitr->halfedge();
    const Point &p0 = he->opposite()->vertex()->point();
    const Point &p1 = he->vertex()->point();
    const Point &p2 = he->next()->vertex()->point();
    const FT area = std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2)));
    const Point barycenter = CGAL::centroid(p0, p1, p2);
    facet_areas.insert(std::pair<Facet_handle, FT>(fitr, area));
    facet_centers.insert(std::pair<Facet_handle, Point>(fitr, barycenter));
  }
  FacetAreaMap area_pmap(facet_areas);
  FacetCenterMap center_pmap(facet_centers);

  // create compact metric approximation algorithm instance
  CompactVSA compact_approx(input,
    get(boost::vertex_point, const_cast<Polyhedron &>(input)));

  // construct metric and fitting functors
  CompactMetric metric(center_pmap);
  PointProxyFitting proxy_fitting(center_pmap, area_pmap);
  compact_approx.set_metric(metric, proxy_fitting);

  // using 200 proxies to approximate the shape
  compact_approx.init_by_number(CompactVSA::Hierarchical, 200);
  for (std::size_t i = 0; i < 30; ++i)
    compact_approx.run_one_step();

  return EXIT_SUCCESS;
}
