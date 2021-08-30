#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Variational_shape_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

typedef Mesh::Property_map<face_descriptor, FT> Face_area_map;
typedef Mesh::Property_map<face_descriptor, Point_3> Face_center_map;
typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;

namespace PMP = CGAL::Polygon_mesh_processing;

// user defined point-wise compact metric
struct Compact_metric_point_proxy {
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
    for(const face_descriptor& f : faces) {
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
  Mesh, Vertex_point_map, Compact_metric_point_proxy> Compact_approx;

/**
 * This file tests the user defined metric.
 */
int main()
{
  Mesh mesh;
  std::ifstream input("./data/sphere.off");
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  const double target_edge_length = 0.05;
  const unsigned int nb_iter = 3;

  std::cout << "Start remeshing. ("
    << std::distance(faces(mesh).first, faces(mesh).second) << " faces)..." << std::endl;
  PMP::isotropic_remeshing(
    faces(mesh),
    target_edge_length,
    mesh,
    PMP::parameters::number_of_iterations(nb_iter));
  std::cout << "Remeshing done. ("
    << std::distance(faces(mesh).first, faces(mesh).second) << " faces)..." << std::endl;

  // construct face normal and area map
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

  // create compact metric approximation algorithm instance
  std::cout << "create compact vas instance" << std::endl;
  Compact_metric_point_proxy error_metric(center_pmap, area_pmap);

  Compact_approx approx(mesh, vpmap, error_metric);

  std::cout << "random seeding and run" << std::endl;
  approx.initialize_seeds(
    CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::RANDOM)
    .max_number_of_proxies(20));
  approx.run(20);
  if (approx.number_of_proxies() != 20)
    return EXIT_FAILURE;

  // extract the approximation
  std::cout << "meshing" << std::endl;
  if (approx.extract_mesh(CGAL::parameters::subdivision_ratio(5.0)))
    std::cout << "manifold." << std::endl;
  else
    std::cout << "non-manifold" << std::endl;

  return EXIT_SUCCESS;
}
