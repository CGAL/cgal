#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Variational_shape_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Point_3 Point;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > Face_area_map;
typedef boost::associative_property_map<std::map<Facet_handle, Point> > Face_center_map;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;

// user defined point-wise compact metric
struct Compact_metric_point_proxy {
// use point as proxy
  typedef Point Proxy;

  // we keep a precomputed property map to speed up computations
  Compact_metric_point_proxy(const Face_center_map &_center_pmap, const Face_area_map &_area_pmap)
    : center_pmap(_center_pmap), area_pmap(_area_pmap) {}

  // compute and return error from a face to a proxy,
  // defined as the Euclidean distance between
  // the face center of mass and proxy point.
  FT compute_error(const Polyhedron &tm, const Facet_handle &f, const Proxy &px) const {
    (void)(tm);
    return FT(std::sqrt(CGAL::to_double(
      CGAL::squared_distance(center_pmap[f], px))));
  }

  // template functor to compute a best-fit 
  // proxy from a range of faces
  template <typename FaceRange>
  Proxy fit_proxy(const FaceRange &faces, const Polyhedron &tm) const {
    (void)(tm);
    // fitting center
    Vector center = CGAL::NULL_VECTOR;
    FT sum_areas = FT(0.0);
    BOOST_FOREACH(const Facet_handle &f, faces) {
      center = center + (center_pmap[f] - CGAL::ORIGIN) * area_pmap[f];
      sum_areas += area_pmap[f];
    }
    center = center / sum_areas;
    return CGAL::ORIGIN + center;
  }

  const Face_center_map center_pmap;
  const Face_area_map area_pmap;
};

typedef CGAL::Variational_shape_approximation<
  Polyhedron, Vertex_point_map, Compact_metric_point_proxy> Compact_approx;

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

  // construct face normal & area map
  std::map<Facet_handle, FT> face_areas;
  std::map<Facet_handle, Point> face_centers;
  for(Facet_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    const Halfedge_handle he = fitr->halfedge();
    const Point &p0 = he->opposite()->vertex()->point();
    const Point &p1 = he->vertex()->point();
    const Point &p2 = he->next()->vertex()->point();
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    face_areas.insert(std::pair<Facet_handle, FT>(fitr, area));
    face_centers.insert(std::pair<Facet_handle, Point>(fitr, CGAL::centroid(p0, p1, p2)));
  }
  Face_area_map area_pmap(face_areas);
  Face_center_map center_pmap(face_centers);

  // create compact metric approximation algorithm instance
  std::cout << "create compact vas instance" << std::endl;
  Compact_metric_point_proxy error_metric(center_pmap, area_pmap);

  Compact_approx approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)),
    error_metric);

  std::cout << "random seeding and run" << std::endl;
  approx.initialize_seeds(CGAL::VSA::parameters::seeding_method(CGAL::VSA::RANDOM)
    .max_nb_of_proxies(20));
  approx.run(20);
  if (approx.proxies_size() != 20)
    return EXIT_FAILURE;

  // extract the approximation
  std::cout << "meshing" << std::endl;
  if (approx.extract_mesh(CGAL::VSA::parameters::subdivision_ratio(5.0)))
    std::cout << "manifold." << std::endl;
  else
    std::cout << "non-manifold" << std::endl;

  return EXIT_SUCCESS;
}
