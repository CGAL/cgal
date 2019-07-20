#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Variational_shape_approximation.h>

namespace VSA = CGAL::Surface_mesh_approximation;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;
typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map> Mesh_approximation;

// L21 error metric
typedef Mesh_approximation::Error_metric L21_metric;

int main()
{
  // reads input surface triangle mesh
  Mesh mesh;
  std::ifstream file("data/bear.off");
  if (!file || !(file >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  Vertex_point_map vpmap = get(boost::vertex_point, const_cast<Mesh &>(mesh));

  // error metric and fitting function
  L21_metric error_metric(mesh, vpmap);

  // creates VSA algorithm instance
  Mesh_approximation approx(mesh, vpmap, error_metric);

  // seeds 100 random proxies
  approx.initialize_seeds(CGAL::parameters::seeding_method(VSA::RANDOM)
      .max_number_of_proxies(100));

  // runs 30 iterations
  approx.run(30);

  // adds 3 proxies to the one with the maximum fitting error,
  // running 5 iterations between each addition
  approx.add_to_furthest_proxies(3, 5);

  // runs 10 iterations
  approx.run(10);

  // teleports 2 proxies to tunnel out of local minima,
  // running 5 iterations between each teleport
  approx.teleport_proxies(2, 5);

  // runs 10 iterations
  approx.run(10);

  // extract approximated mesh with default parameters
  approx.extract_mesh(CGAL::parameters::all_default());

  // get approximated triangle soup
  std::vector<Kernel::Point_3> anchors;
  std::vector<std::array<std::size_t, 3> > triangles;
  approx.output(CGAL::parameters::anchors(std::back_inserter(anchors)).
    triangles(std::back_inserter(triangles)));

  return EXIT_SUCCESS;
}
