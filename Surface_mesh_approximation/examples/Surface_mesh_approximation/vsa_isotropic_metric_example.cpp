#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

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

typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > Facet_area_map;
typedef boost::associative_property_map<std::map<Facet_handle, Point> > Facet_center_map;

// user-defined "compact" error metric using type Point_3 as proxy
struct Compact_metric_point_proxy
{
  // use point as proxy
  typedef Point Proxy;

  // we keep a precomputed property map to speed up computations
  Compact_metric_point_proxy(const Facet_center_map &_center_pmap, const Facet_area_map &_area_pmap)
    : center_pmap(_center_pmap), area_pmap(_area_pmap) {}

  // compute and return error from a facet to a proxy,
  // defined as the Euclidean distance between
  // the facet center of mass and proxy point.
  FT compute_error(const Polyhedron &tm, const Facet_handle &f, const Proxy &px) const {
    (void)(tm);
    return FT(std::sqrt(CGAL::to_double(
      CGAL::squared_distance(center_pmap[f], px))));
  }

  // template functor to compute a best-fit 
  // proxy from a range of faces
  template <typename FaceRange>
  Proxy fit_proxy(const Polyhedron &tm, const FaceRange &faces) const {
    (void)(tm);
    // fitting center
    Vector center = CGAL::NULL_VECTOR;
    FT sum_areas = FT(0.0);
    BOOST_FOREACH(const Facet_handle &f, faces) {
      center = center + (center_pmap[f] - CGAL::ORIGIN) * area_pmap[f];
      sum_areas += area_pmap[f];
    }
    center = center / sum_areas; // TODO: deal with case where sum = 0
    return CGAL::ORIGIN + center;
  }

  const Facet_center_map center_pmap;
  const Facet_area_map area_pmap;
};

typedef CGAL::VSA_approximation<
  Polyhedron, Vertex_point_map, Compact_metric_point_proxy> Approximation;

int main()
{
  // reads input mesh
  Polyhedron input;
  std::ifstream file("data/bear.off");
  file >> input;

  // constructs precomputed facet normal and area map
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
  Facet_area_map area_pmap(facet_areas);
  Facet_center_map center_pmap(facet_centers);

  // creates compact metric approximation algorithm instance
  Approximation approx(input,
    get(boost::vertex_point, const_cast<Polyhedron &>(input)));

  // constructs and set metric
  Compact_metric_point_proxy metric(center_pmap, area_pmap);
  approx.set_metric(metric);

  // approximates via 200 proxies and 30 iterations
  approx.seeding(CGAL::Hierarchical, 200);
  approx.run(30);

  return EXIT_SUCCESS;
}
