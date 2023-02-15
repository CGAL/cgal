#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/Surface_mesh_approximation/L2_metric_plane_proxy.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

typedef Mesh::Property_map<face_descriptor, std::size_t> Face_proxy_map;
typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::Surface_mesh_approximation::L2_metric_plane_proxy<Mesh> L2_metric_plane_proxy;
typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map, L2_metric_plane_proxy> L2_approx;
typedef L2_approx::Proxy Plane_proxy;

namespace PMP = CGAL::Polygon_mesh_processing;

/**
 * This file tests the main class API and the L2 metric.
 * It should cover all the APIs.
 */
int main()
{
  Mesh mesh;
  std::ifstream input(CGAL::data_file_path("meshes/sphere.off"));
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
    CGAL::parameters::number_of_iterations(nb_iter));
  std::cout << "Remeshing done. ("
    << std::distance(faces(mesh).first, faces(mesh).second) << " faces)..." << std::endl;

  // face proxy index map
  Face_proxy_map fpxmap =
    mesh.add_property_map<face_descriptor, std::size_t>("f:proxy_id", 0).first;

  // create L2_approx L2 metric approximation algorithm instance
  std::cout << "setup algorithm instance" << std::endl;
  L2_metric_plane_proxy error_metric(mesh,
    get(boost::vertex_point, const_cast<Mesh &>(mesh)));
  L2_approx approx(mesh,
    get(boost::vertex_point, const_cast<Mesh &>(mesh)),
    error_metric);

  // random seeding and run
  std::cout << "random seeding and run" << std::endl;
  approx.initialize_seeds(CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::RANDOM)
    .max_number_of_proxies(10));
  approx.run(10);
  if (approx.number_of_proxies() != 10)
    return EXIT_FAILURE;

  // incremental add and run to convergence
  std::cout << "incremental add and run to convergence" << std::endl;
  approx.add_to_furthest_proxies(3, 5);
  if (approx.run_to_convergence(0.1))
    std::cout << "Converged." << std::endl;
  if (approx.number_of_proxies() != 13)
    return EXIT_FAILURE;

  std::cout << "hierarchical add and run" << std::endl;
  approx.add_proxies_error_diffusion(3);
  approx.run(10);
  if (approx.number_of_proxies() != 16)
    return EXIT_FAILURE;

  // merge and teleport the proxies from local minimal
  std::cout << "teleport" << std::endl;
  approx.teleport_proxies(3);
  approx.run(10);
  if (approx.number_of_proxies() != 16)
    return EXIT_FAILURE;

  // split proxy 0 into 2 proxies
  // precondition: proxy 0 should have more than 2 faces
  std::cout << "splitting" << std::endl;
  if (!approx.split(0, 2, 10))
    return EXIT_FAILURE;
  if (approx.number_of_proxies() != 17)
    return EXIT_FAILURE;

  // extract the approximation Mesh
  std::cout << "meshing" << std::endl;
  if (approx.extract_mesh(CGAL::parameters::subdivision_ratio(1.0)))
    std::cout << "manifold." << std::endl;
  else
    std::cout << "non-manifold" << std::endl;

  // get outputs
  std::cout << "get outputs" << std::endl;
  approx.proxy_map(fpxmap);

  for (std::size_t i = 0; i < approx.number_of_proxies(); ++i) {
    std::list<face_descriptor> patch;
    approx.proxy_region(i, std::back_inserter(patch));
  }

  std::vector<Plane_proxy> proxies;
  approx.proxies(std::back_inserter(proxies));

  std::vector<Kernel::Point_3> anchor_pos;
  approx.anchor_points(std::back_inserter(anchor_pos));

  std::vector<vertex_descriptor> anchor_vtx;
  approx.anchor_vertices(std::back_inserter(anchor_vtx));

  std::vector<std::array<std::size_t, 3> > tris;
  approx.indexed_triangles(std::back_inserter(tris));

  std::vector<std::vector<std::size_t> > boundary;
  approx.indexed_boundary_polygons(std::back_inserter(boundary));

  const Kernel::FT drop(0.001);
  const std::size_t iterations = 5;
  std::cout << "re-initialize and hierarchical seeding" << std::endl;
  approx.initialize_seeds(CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::HIERARCHICAL)
    .min_error_drop(drop)
    .number_of_relaxations(iterations));
  approx.run(10);
  std::cout << "#proxies " << approx.number_of_proxies() << std::endl;

  std::cout << "re-initialize and incremental seeding" << std::endl;
  approx.initialize_seeds(CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::INCREMENTAL)
    .min_error_drop(drop)
    .number_of_relaxations(iterations));
  approx.run(10);
  std::cout << "#proxies " << approx.number_of_proxies() << std::endl;

  // extract the approximation Mesh
  std::cout << "meshing" << std::endl;
  if (approx.extract_mesh(CGAL::parameters::subdivision_ratio(1.0)))
    std::cout << "manifold." << std::endl;
  else
    std::cout << "non-manifold" << std::endl;

  // get outputs
  std::cout << "get outputs" << std::endl;
  proxies.clear();
  anchor_pos.clear();
  tris.clear();
  approx.output(CGAL::parameters::face_proxy_map(fpxmap).
    proxies(std::back_inserter(proxies)).
    anchors(std::back_inserter(anchor_pos)).
    triangles(std::back_inserter(tris)));

  return EXIT_SUCCESS;
}
