#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Simple_cartesian<double> Sckernel;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename TM>
int load_and_remesh_sm(TM &mesh) {
  std::ifstream input("./data/cube.off");
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

  return EXIT_SUCCESS;
}

template <typename TM>
int load_and_remesh_poly(TM &mesh) {
  std::ifstream input("./data/cube.off");
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
    PMP::parameters::number_of_iterations(nb_iter).
    face_index_map(get(boost::face_external_index, mesh)));
  std::cout << "Remeshing done. ("
    << std::distance(faces(mesh).first, faces(mesh).second) << " faces)..." << std::endl;

  return EXIT_SUCCESS;
}

template <typename K, typename TM>
void run_approximation(const TM &mesh) {
  std::vector<typename K::Point_3> points;
  std::vector<std::array<std::size_t, 3> > triangles;

  CGAL::Surface_mesh_approximation::approximate_triangle_mesh(mesh,
    CGAL::parameters::max_number_of_proxies(6).
      number_of_iterations(30).
      number_of_relaxations(5).
      subdivision_ratio(0.5).
      anchors(std::back_inserter(points)).
      triangles(std::back_inserter(triangles)));
}

template <typename K, typename TM>
int test_sm() {
  TM mesh;
  if (load_and_remesh_sm(mesh) == EXIT_FAILURE)
    return EXIT_FAILURE;
  run_approximation<K, TM>(mesh);
  return EXIT_SUCCESS;
}

template <typename K, typename TM>
int test_poly() {
  TM mesh;
  if (load_and_remesh_poly(mesh) == EXIT_FAILURE)
    return EXIT_FAILURE;
  run_approximation<K, TM>(mesh);
  return EXIT_SUCCESS;
}

/**
 * This file tests the VSA free function API with different configuration of kernel and surface mesh.
 */
int main()
{
  if (test_poly<Epic, CGAL::Polyhedron_3<Epic> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (test_sm<Epic, CGAL::Surface_mesh<Epic::Point_3> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (test_poly<Sckernel, CGAL::Polyhedron_3<Sckernel> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (test_sm<Sckernel, CGAL::Surface_mesh<Sckernel::Point_3> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
