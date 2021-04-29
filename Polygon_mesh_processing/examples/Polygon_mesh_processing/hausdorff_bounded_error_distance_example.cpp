#include <CGAL/Bbox_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/PLY.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

using TAG                     = CGAL::Parallel_if_available_tag;
using Surface_mesh            = CGAL::Surface_mesh<Point_3>;
using Polyhedron              = CGAL::Polyhedron_3<Kernel>;
using Affine_transformation_3 = CGAL::Aff_transformation_3<Kernel>;
using Timer                   = CGAL::Real_timer;

namespace PMP = CGAL::Polygon_mesh_processing;

template<class Mesh>
void get_mesh(const std::string filepath, Mesh& mesh) {

  mesh.clear();
  std::ifstream input(filepath);
  input >> mesh;
  std::cout << "* getting mesh with " << num_faces(mesh) << " faces" << std::endl;
}

void get_meshes(
  const std::string filepath1, const std::string filepath2,
  Surface_mesh& mesh1, Surface_mesh& mesh2) {

  get_mesh(filepath1, mesh1);
  get_mesh(filepath2, mesh2);
}

void save_mesh(const Surface_mesh& mesh, const std::string filepath) {

  if (!CGAL::write_PLY(filepath + ".ply", mesh)) {
    std::cerr << "ERROR: cannot save this file: " << filepath << std::endl;
    exit(EXIT_FAILURE);
  }
}

// An easy example of a tetrahedron and its remeshed version.
void remeshing_tetrahedon_example(
  const double error_bound, const bool save = false) {

  Timer timer;
  Surface_mesh mesh1, mesh2;
  std::cout << std::endl << "* (E1) remeshing tetrahedron example:" << std::endl;

  CGAL::make_tetrahedron(
    Point_3(0, 0, 0), Point_3(2, 0, 0),
    Point_3(1, 1, 1), Point_3(1, 0, 2), mesh1);
  mesh2 = mesh1;

  using edge_descriptor = typename boost::graph_traits<Surface_mesh>::edge_descriptor;
  Surface_mesh::Property_map<edge_descriptor, bool> is_constrained_map =
    mesh2.add_property_map<edge_descriptor, bool>("e:is_constrained", true).first;

  const double target_edge_length = 0.05;
  PMP::isotropic_remeshing(
    mesh2.faces(), target_edge_length, mesh2,
    PMP::parameters::edge_is_constrained_map(is_constrained_map));

  if (save) save_mesh(mesh1, "mesh1");
  if (save) save_mesh(mesh2, "mesh2");

  timer.reset();
  timer.start();
  std::cout << "* bounded Hausdorff distance: " <<
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound) << std::endl;
  timer.stop();
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;
}

// Example with a point realizing the Hausdorff distance strictly
// lying in the interior of a triangle.
void interior_triangle_example(
  const double error_bound, const bool save = false) {

  Timer timer;
  Surface_mesh mesh1, mesh2;
  std::cout << std::endl << "* (E2) interior triangle example:" << std::endl;

  mesh1.add_vertex(Point_3(-1,  1, 1));
  mesh1.add_vertex(Point_3( 0, -1, 1));
  mesh1.add_vertex(Point_3( 1,  1, 1));
  mesh1.add_face(mesh1.vertices());

  auto v0 = mesh2.add_vertex(Point_3(-1.0,  1,  0));
  auto v1 = mesh2.add_vertex(Point_3( 0.0, -1,  0));
  auto v2 = mesh2.add_vertex(Point_3( 1.0,  1,  0));
  auto v3 = mesh2.add_vertex(Point_3( 0.0,  1, -1));
  auto v4 = mesh2.add_vertex(Point_3(-0.5,  0, -1));
  auto v5 = mesh2.add_vertex(Point_3( 0.5,  0, -1));
  mesh2.add_face(v0, v3, v4);
  mesh2.add_face(v1, v4, v5);
  mesh2.add_face(v2, v5, v3);

  if (save) save_mesh(mesh1, "mesh1");
  if (save) save_mesh(mesh2, "mesh2");

  timer.reset();
  timer.start();
  std::cout << "* bounded Hausdorff distance: " <<
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound) << std::endl;
  timer.stop();
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;
}

// Read a real mesh given by the user, perturb it slightly, and compute the
// Hausdorff distance between the original mesh and its pertubation.
void perturbing_mesh_example(
  const std::string filepath, const double error_bound, const bool save = false) {

  Timer timer;
  std::cout << std::endl << "* (E3) perturbing mesh example:" << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath, filepath, mesh1, mesh2);

  const double max_size = 0.1;
  PMP::random_perturbation(
    mesh2.vertices(), mesh2, max_size, CGAL::parameters::do_project(false));
  std::cout << "* perturbing the second mesh" << std::endl;

  if (save) save_mesh(mesh2, "mesh2");

  timer.reset();
  timer.start();
  std::cout << "* bounded Hausdorff distance: " <<
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound) << std::endl;
  timer.stop();
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;
}

// Read two meshes given by the user, initially place them at their originally
// given position. Move the second mesh in 300 steps away from the first one.
// Print how the Hausdorff distance changes.
void moving_mesh_example(
  const std::string filepath1, const std::string filepath2,
  const std::size_t n, const double error_bound, const bool save = false) {

  Timer timer;
  std::cout << std::endl << "* (E4) moving mesh example:" << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath1, filepath2, mesh1, mesh2);

  const auto bbox = PMP::bbox(mesh2);
  const FT distance = static_cast<FT>(CGAL::sqrt(CGAL::to_double(
    CGAL::squared_distance(
      Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin()),
      Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax())))));

  const FT t = FT(1) / FT(100);
  if (save) save_mesh(mesh2, "mesh-0");

  for (std::size_t i = 0; i < n; ++i) {
    PMP::transform(Affine_transformation_3(CGAL::Translation(),
      Vector_3(t * distance, t * distance, t * distance)), mesh2);

    timer.reset();
    timer.start();
    std::cout << "* position: " << i << std::endl;
    std::cout << "* bounded Hausdorff distance: " <<
      PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound) << std::endl;
    timer.stop();
    std::cout << "* processing time: " << timer.time() << " s." << std::endl;
    if (save) save_mesh(mesh2, "mesh-" + std::to_string(i + 1));
  }
}

void perturbing_mesh_example_with_polyhedron(
  const std::string filepath, const double error_bound) {

  Timer timer;
  std::cout << std::endl << "* (E3) perturbing mesh example:" << std::endl;

  Surface_mesh mesh1;
  Polyhedron mesh2;
  get_mesh(filepath, mesh1);
  get_mesh(filepath, mesh2);

  const double max_size = 0.1;
  PMP::random_perturbation(
    vertices(mesh2), mesh2, max_size, CGAL::parameters::do_project(false));
  std::cout << "* perturbing the second mesh" << std::endl;

  timer.reset();
  timer.start();
  std::cout << "* bounded Hausdorff distance 1->2: " <<
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound) << std::endl;
  std::cout << "* bounded Hausdorff distance 2->1: " <<
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh2, mesh1, error_bound) << std::endl;
  timer.stop();
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;
}

int main(int argc, char** argv) {

  const double error_bound = 1e-4;
  std::cout << std::endl << "* error bound: " << error_bound << std::endl;
  std::string filepath = (argc > 1 ? argv[1] : "data/blobby.off");

  remeshing_tetrahedon_example(error_bound);
  interior_triangle_example(error_bound);
  perturbing_mesh_example(filepath, error_bound);
  perturbing_mesh_example_with_polyhedron(filepath, error_bound);
  moving_mesh_example(filepath, filepath, 5, error_bound);

  std::cout << std::endl;
  return EXIT_SUCCESS;
}
