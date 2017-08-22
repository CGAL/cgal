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
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Polyhedron_3::Facet_handle Facet_handle;
typedef Polyhedron_3::Halfedge_handle Halfedge_handle;
typedef Polyhedron_3::Facet_iterator Facet_iterator;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;
typedef boost::associative_property_map<std::map<Facet_handle, Point_3> > FacetCenterMap;

// proxy
struct PointProxy {
  Facet_handle seed;
  Point_3 center;
};

// metric functor
struct CompactMetric {
  typedef PointProxy Proxy;

  // keeping a precomputed property map to speed up the computation
  CompactMetric(const FacetCenterMap &_center_pmap)
    : center_pmap(_center_pmap) {}

  // calculate the error of a facet to a proxy
  // here is just the Euclidean distance
  FT operator()(const Facet_handle &f, const PointProxy &px) const {
    return FT(std::sqrt(CGAL::to_double(
      CGAL::squared_distance(center_pmap[f], px.center))));
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
  PointProxy operator()(const FacetIterator beg, const FacetIterator end) const {
    // fitting center
    Vector_3 center = CGAL::NULL_VECTOR;
    FT area(0);
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      center = center + (center_pmap[*fitr] - CGAL::ORIGIN) * area_pmap[*fitr];
      area += area_pmap[*fitr];
    }
    center = center / area;

    // construct proxy
    PointProxy px;
    px.center = CGAL::ORIGIN + center;

    return px;
  }

  const FacetCenterMap center_pmap;
  const FacetAreaMap area_pmap;
};
typedef CGAL::VSA_approximation<Polyhedron_3, PointProxy, CompactMetric, PointProxyFitting> CompactVSA;

int main()
{
  // create and read Polyhedron_3
  Polyhedron_3 mesh;
  std::ifstream input("data/bear.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // construct precomputed facet normal and area map
  std::map<Facet_handle, FT> facet_areas;
  std::map<Facet_handle, Point_3> facet_centers;
  for(Facet_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    const Halfedge_handle he = fitr->halfedge();
    const Point_3 &p0 = he->opposite()->vertex()->point();
    const Point_3 &p1 = he->vertex()->point();
    const Point_3 &p2 = he->next()->vertex()->point();
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<Facet_handle, FT>(fitr, area));
    facet_centers.insert(std::pair<Facet_handle, Point_3>(fitr, CGAL::centroid(p0, p1, p2)));
  }
  FacetAreaMap area_pmap(facet_areas);
  FacetCenterMap center_pmap(facet_centers);

  CompactMetric metric(center_pmap);
  PointProxyFitting proxy_fitting(center_pmap, area_pmap);

  // create compact metric approximation algorithm instance
  CompactVSA compact_approx;
  compact_approx.set_mesh(mesh);
  compact_approx.set_error_metric(metric);
  compact_approx.set_proxy_fitting(proxy_fitting);

  // using 200 proxies to approximate the shape
  compact_approx.init_proxies(200, CompactVSA::HierarchicalInit);
  for (std::size_t i = 0; i < 30; ++i)
    compact_approx.run_one_step();

  return EXIT_SUCCESS;
}
