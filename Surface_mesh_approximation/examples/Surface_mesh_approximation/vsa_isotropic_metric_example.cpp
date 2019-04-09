#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/property_map.h>
#include <CGAL/Variational_shape_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;
typedef Mesh::Property_map<face_descriptor, FT> Face_area_map;
typedef Mesh::Property_map<face_descriptor, Point_3> Face_center_map;

namespace VSA = CGAL::Surface_mesh_approximation;

// user-defined "compact" error metric using type Point_3 as proxy
struct Compact_metric_point_proxy
{
  // use point as proxy
  typedef Point_3 Proxy;

  // we keep a precomputed property map to speed up computations
  Compact_metric_point_proxy(const Face_center_map &center_pmap_, const Face_area_map &area_pmap_)
    : center_pmap(center_pmap_), area_pmap(area_pmap_) {}

  // compute and return error from a face to a proxy,
  // defined as the Euclidean distance between
  // the face center of mass and proxy point.
  FT compute_error(const face_descriptor &f, const Mesh &, const Proxy &px) const {
    return CGAL::sqrt(CGAL::squared_distance(center_pmap[f], px));
  }

  // template functor to compute a best-fit
  // proxy from a range of faces
  template <typename FaceRange>
  Proxy fit_proxy(const FaceRange &faces, const Mesh &) const {
    // fitting center
    Vector_3 center = CGAL::NULL_VECTOR;
    FT sum_areas = FT(0.0);
    for(const face_descriptor f : faces) {
      center = center + (center_pmap[f] - CGAL::ORIGIN) * area_pmap[f];
      sum_areas += area_pmap[f];
    }
    // deal with case where sum = 0
    if (center == CGAL::NULL_VECTOR || sum_areas <= FT(0.0)) {
      std::cerr << "Error: degenerate geometry." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    center = center / sum_areas;
    return CGAL::ORIGIN + center;
  }

  const Face_center_map center_pmap;
  const Face_area_map area_pmap;
};

typedef CGAL::Variational_shape_approximation<
  Mesh, Vertex_point_map, Compact_metric_point_proxy> Approximation;

int main()
{
  // reads input mesh
  Mesh mesh;
  std::ifstream file("data/bear.off");
  file >> mesh;

  // constructs precomputed face normal and area map
  Vertex_point_map vpmap = get(boost::vertex_point, const_cast<Mesh &>(mesh));
  Face_area_map area_pmap =
    mesh.add_property_map<face_descriptor, FT>("f:area", FT(0.0)).first;
  Face_center_map center_pmap =
    mesh.add_property_map<face_descriptor, Point_3>("f:center", CGAL::ORIGIN).first;
  for(face_descriptor f : faces(mesh)) {
    const halfedge_descriptor he = halfedge(f, mesh);
    const Point_3 &p0 = vpmap[source(he, mesh)];
    const Point_3 &p1 = vpmap[target(he, mesh)];
    const Point_3 &p2 = vpmap[target(next(he, mesh), mesh)];
    put(area_pmap, f, CGAL::sqrt(CGAL::squared_area(p0, p1, p2)));
    put(center_pmap, f, CGAL::centroid(p0, p1, p2));
  }

  // error metric and fitting function
  Compact_metric_point_proxy error_metric(center_pmap, area_pmap);

  // creates compact metric approximation algorithm instance
  Approximation approx(mesh, vpmap, error_metric);

  // approximates via 200 proxies and 30 iterations
  approx.initialize_seeds(
    CGAL::parameters::seeding_method(VSA::HIERARCHICAL)
    .max_number_of_proxies(200));
  approx.run(30);

  return EXIT_SUCCESS;
}
