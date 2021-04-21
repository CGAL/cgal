#include <CGAL/Bbox_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/PLY.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

using SCK   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel   = EPICK;
using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

using TAG                     = CGAL::Parallel_if_available_tag;
using Surface_mesh            = CGAL::Surface_mesh<Point_3>;
using Affine_transformation_3 = CGAL::Aff_transformation_3<Kernel>;
using Timer                   = CGAL::Real_timer;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Approximate_hd_wrapper {
  const std::size_t m_num_samples;
  std::string name() const { return "approximate"; }
  Approximate_hd_wrapper(const std::size_t num_samples) : m_num_samples(num_samples) { }
  double operator()(const Surface_mesh& tm1, const Surface_mesh& tm2) const {
    return PMP::approximate_Hausdorff_distance<TAG>(tm1, tm2,
      PMP::parameters::number_of_points_per_area_unit(m_num_samples),
      PMP::parameters::number_of_points_per_area_unit(m_num_samples));
  }
};

struct Naive_bounded_error_hd_wrapper {
  const double m_error_bound;
  std::string name() const { return "naive bounded error"; }
  Naive_bounded_error_hd_wrapper(const double error_bound) : m_error_bound(error_bound) { }
  double operator()(const Surface_mesh& tm1, const Surface_mesh& tm2) const {
    return PMP::bounded_error_Hausdorff_distance_naive<TAG>(tm1, tm2, m_error_bound);
  }
};

struct Bounded_error_hd_wrapper {
  const double m_error_bound;
  std::string name() const { return "bounded error"; }
  Bounded_error_hd_wrapper(const double error_bound) : m_error_bound(error_bound) { }
  double operator()(const Surface_mesh& tm1, const Surface_mesh& tm2) const {
    return PMP::bounded_error_Hausdorff_distance<TAG>(tm1, tm2, m_error_bound);
  }
};

void get_mesh(const std::string filepath, Surface_mesh& mesh) {

  mesh.clear();
  std::ifstream input(filepath);
  input >> mesh;
  std::cout << "* getting mesh with " << mesh.number_of_faces() << " faces" << std::endl;
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
void remeshing_tetrahedon_example(const double error_bound) {

  Timer timer;
  Surface_mesh mesh1, mesh2;
  std::cout << std::endl << "* (E1) remeshing tetrahedron example:" << std::endl;

  CGAL::make_tetrahedron(
    Point_3(0, 0, 0), Point_3(2, 0, 0),
    Point_3(1, 1, 1), Point_3(1, 0, 2), mesh1);
  mesh2 = mesh1;

  // TODO: How to preserve edges?
  const double target_edge_length = 0.05;
  PMP::isotropic_remeshing(mesh2.faces(), target_edge_length, mesh2);

  // save_mesh(mesh1, "1-mesh1");
  // save_mesh(mesh2, "1-mesh2");

  timer.reset();
  timer.start();
  std::cout << "* bounded Hausdorff distance: " <<
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound) << std::endl;
  timer.stop();
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;
}

// Example with a point realizing the Hausdorff distance strictly
// lying in the interior of a triangle.
void interior_triangle_example(const double error_bound) {

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

  // save_mesh(mesh1, "2-mesh1");
  // save_mesh(mesh2, "2-mesh2");

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
  const std::string filepath, const double error_bound) {

  Timer timer;
  std::cout << std::endl << "* (E3) perturbing mesh example:" << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath, filepath, mesh1, mesh2);

  const double max_size = 0.1;
  PMP::random_perturbation(
    mesh2.vertices(), mesh2, max_size, CGAL::parameters::do_project(false));
  std::cout << "* perturbing the second mesh" << std::endl;

  // save_mesh(mesh2, "3-mesh1");
  // save_mesh(mesh2, "3-mesh2");

  timer.reset();
  timer.start();
  std::cout << "* bounded Hausdorff distance: " <<
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh2, mesh2, error_bound) << std::endl;
  timer.stop();
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;
}

// Read two meshes given by the user, initially place them at their originally
// given position. Move the second mesh in 300 steps away from the first one.
// Print how the Hausdorff distance changes.
void moving_mesh_example(
  const std::string filepath1, const std::string filepath2,
  const std::size_t n, const double error_bound) {

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
  // save_mesh(mesh2, "4-mesh-0");

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
    // save_mesh(mesh2, "4-mesh-" + std::to_string(i + 1));
  }
}

template<typename FunctionWrapper>
void test_0(const FunctionWrapper& functor, const bool save = false) {

  // The same triangle.
  // Expected distance is 0.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 0 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  mesh2 = mesh1;
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 0): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 0): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_1(const FunctionWrapper& functor, const bool save = false) {

  // Two triangles are parallel and 1 unit distance away from each other.
  // Expected distance is 1.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 1 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  mesh2.add_vertex(Point_3(0, 0, 1));
  mesh2.add_vertex(Point_3(2, 0, 1));
  mesh2.add_vertex(Point_3(1, 1, 1));
  mesh2.add_face(mesh2.vertices());
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 1): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 1): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_2(const FunctionWrapper& functor, const bool save = false) {

  // One triangle is orthogonal to the other one and shares a common edge.
  // Expected distance is 1.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 2 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  mesh2.add_vertex(Point_3(0, 0, 0));
  mesh2.add_vertex(Point_3(2, 0, 0));
  mesh2.add_vertex(Point_3(1, 0, 1));
  mesh2.add_face(mesh2.vertices());
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 1): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 1): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_3(const FunctionWrapper& functor, const bool save = false) {

  // One triangle is orthogonal to the other one and shares a common edge
  // that is moved 1 unit distance away.
  // Expected distances are sqrt(2) and 2.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 3 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  mesh2.add_vertex(Point_3(0, 0, 1));
  mesh2.add_vertex(Point_3(2, 0, 1));
  mesh2.add_vertex(Point_3(1, 0, 2));
  mesh2.add_face(mesh2.vertices());
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected sqrt(2)): " << dista << std::endl;
  std::cout << "* HInverted distance (expected      2 ): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_4(const FunctionWrapper& functor, const bool save = false) {

  // One triangle is orthogonal to the other one and shares a common vertex.
  // Expected distance is 1.2247448713915889407.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 4 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  mesh2.add_vertex(Point_3(0, 1, 1));
  mesh2.add_vertex(Point_3(2, 1, 1));
  mesh2.add_vertex(Point_3(1, 1, 0));
  mesh2.add_face(mesh2.vertices());
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 1.22): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 1.22): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_5(const FunctionWrapper& functor, const bool save = false) {

  // One triangle is orthogonal to the other one and shares a common vertex
  // that is moved 1 unit distance away.
  // Expected distances are 1.7320508075688771932 and 2.1213203435596423851.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 5 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  mesh2.add_vertex(Point_3(0, 1, 2));
  mesh2.add_vertex(Point_3(2, 1, 2));
  mesh2.add_vertex(Point_3(1, 1, 1));
  mesh2.add_face(mesh2.vertices());
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 1.73): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 2.12): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_6(const FunctionWrapper& functor, const bool save = false) {

  // The first and second mesh have different number of triangles.
  // They are parallel and lie at the same plane. The middle triangle is overlapping.
  // Expected distances are 0 and 0.70710678118654757274.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 6 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  const auto v0 = mesh2.add_vertex(Point_3(0, 0, 0));
  const auto v1 = mesh2.add_vertex(Point_3(2, 0, 0));
  const auto v2 = mesh2.add_vertex(Point_3(1, 1, 0));
  const auto v3 = mesh2.add_vertex(Point_3(2, 1, 0));
  const auto v4 = mesh2.add_vertex(Point_3(0, 1, 0));

  mesh2.add_face(v0, v1, v2);
  mesh2.add_face(v2, v1, v3);
  mesh2.add_face(v0, v2, v4);
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 0.0): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 0.7): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_7(const FunctionWrapper& functor, const bool save = false) {

  // One triangle is moved to 0.5 unit distance away from 3 other triangles.
  // The first and second meshes are parallel.
  // Expected distances are 0.5 and 0.86602540378443859659.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 7 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  const auto v0 = mesh2.add_vertex(Point_3(0, 0, 0.5));
  const auto v1 = mesh2.add_vertex(Point_3(2, 0, 0.5));
  const auto v2 = mesh2.add_vertex(Point_3(1, 1, 0.5));
  const auto v3 = mesh2.add_vertex(Point_3(2, 1, 0.5));
  const auto v4 = mesh2.add_vertex(Point_3(0, 1, 0.5));

  mesh2.add_face(v0, v1, v2);
  mesh2.add_face(v2, v1, v3);
  mesh2.add_face(v0, v2, v4);
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 0.50): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 0.86): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_8(const FunctionWrapper& functor, const bool save = false) {

  // One mesh has one triangle at zero level, another mesh has two triangles
  // where the first one is at level 1 and the second one is at level 2.
  // Expected distances are 1 and 2.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 8 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  if (save) save_mesh(mesh1, "mesh1");

  const auto v0 = mesh2.add_vertex(Point_3(0, 0, 1));
  const auto v1 = mesh2.add_vertex(Point_3(2, 0, 1));
  const auto v2 = mesh2.add_vertex(Point_3(1, 1, 1));
  const auto v3 = mesh2.add_vertex(Point_3(0, 0, 2));
  const auto v4 = mesh2.add_vertex(Point_3(2, 0, 2));
  const auto v5 = mesh2.add_vertex(Point_3(1, 1, 2));

  mesh2.add_face(v0, v1, v2);
  mesh2.add_face(v3, v4, v5);
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 1): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 2): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_9(const FunctionWrapper& functor, const bool save = false) {

  // Two meshes partially overlap, have 2 triangles in common and each one has
  // two its own trianles. All triangles form a Z shape where the height is 1.
  // The expected result is 1.

  std::cout.precision(20);
  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 9 ---- " << std::endl;

  auto v0 = mesh1.add_vertex(Point_3(0, 0, 0));
  auto v1 = mesh1.add_vertex(Point_3(1, 0, 0));
  auto v2 = mesh1.add_vertex(Point_3(0, 1, 0));
  auto v3 = mesh1.add_vertex(Point_3(1, 1, 0));
  auto v4 = mesh1.add_vertex(Point_3(1, 0, 1));
  auto v5 = mesh1.add_vertex(Point_3(1, 1, 1));
  mesh1.add_face(v0, v1, v2);
  mesh1.add_face(v2, v1, v3);
  mesh1.add_face(v1, v4, v3);
  mesh1.add_face(v3, v4, v5);
  if (save) save_mesh(mesh1, "mesh1");

  v0 = mesh2.add_vertex(Point_3(2, 0, 1));
  v1 = mesh2.add_vertex(Point_3(1, 0, 0));
  v2 = mesh2.add_vertex(Point_3(2, 1, 1));
  v3 = mesh2.add_vertex(Point_3(1, 1, 0));
  v4 = mesh2.add_vertex(Point_3(1, 0, 1));
  v5 = mesh2.add_vertex(Point_3(1, 1, 1));
  mesh2.add_face(v1, v4, v3);
  mesh2.add_face(v3, v4, v5);
  mesh2.add_face(v4, v0, v5);
  mesh2.add_face(v5, v0, v2);
  if (save) save_mesh(mesh2, "mesh2");

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);

  std::cout << "* Hausdorff distance (expected 1): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 1): " << distb << std::endl;
}

template<typename FunctionWrapper>
void test_synthetic_data(const FunctionWrapper& functor) {

  std::cout << std::endl << "* testing synthetic data:" << std::endl;
  std::cout << "* name -> " << functor.name() << std::endl;

  test_0(functor); // 1 parallel
  test_1(functor);
  test_2(functor); // 1 edge touching
  test_3(functor);
  test_4(functor); // 1 vertex touching
  test_5(functor);
  test_6(functor); // 1 to multiple
  test_7(functor);
  test_8(functor);
  test_9(functor); // overlapping
}

template<
typename FunctionWrapper1,
typename FunctionWrapper2>
void test_one_versus_another(
  const FunctionWrapper1& functor1,
  const FunctionWrapper2& functor2) {

  std::cout.precision(20);
  const std::string filepath1 = "data/tetrahedron.off";
  const std::string filepath2 = "data/tetrahedron-remeshed.off";

  std::cout << std::endl << "* testing one versus another (tetrahedron):" << std::endl;
  std::cout << "* name 1 -> " << functor1.name() << std::endl;
  std::cout << "* name 2 -> " << functor2.name() << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath1, filepath2, mesh1, mesh2);

  // TEST 0.
  // Load and compare.
  // The expected distance is 0.
  std::cout << " ---- testing 0 ---- " << std::endl;
  const double dista0 = functor1(mesh1, mesh2);
  const double distb0 = functor2(mesh1, mesh2);
  std::cout << "* Hausdorff distance1: " << dista0 << std::endl;
  std::cout << "* Hausdorff distance2: " << distb0 << std::endl;

  std::cout << "* traslating by 1 unit ..." << std::endl;
  PMP::transform(Affine_transformation_3(CGAL::Translation(),
    Vector_3(FT(0), FT(0), FT(1))), mesh2);

  // TEST 1.
  // Translate by 1 unit distance and compare.
  // The expected distance is 1.
  std::cout << " ---- testing 1 ---- " << std::endl;
  const double dista1 = functor1(mesh1, mesh2);
  const double distb1 = functor2(mesh1, mesh2);
  std::cout << "* Hausdorff distance1: " << dista1 << std::endl;
  std::cout << "* Hausdorff distance2: " << distb1 << std::endl;
}

template<
typename FunctionWrapper1,
typename FunctionWrapper2>
void test_real_meshes(
  const std::string filepath1,
  const std::string filepath2,
  const FunctionWrapper1& functor1,
  const FunctionWrapper2& functor2) {

  std::cout.precision(20);
  std::cout << std::endl << "* testing real meshes:" << std::endl;
  std::cout << "* input path 1: " << filepath1 << std::endl;
  std::cout << "* input path 2: " << filepath2 << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath1, filepath2, mesh1, mesh2);

  std::cout << std::endl;
  std::cout << "* name 1 -> " << functor1.name() << std::endl;
  std::cout << "* name 2 -> " << functor2.name() << std::endl;

  // Load and compare.
  std::cout << std::endl;
  std::cout << " ---- testing ---- " << std::endl;
  const double dista0 = functor1(mesh1, mesh2);
  const double dista1 = functor1(mesh2, mesh1);
  const double distb0 = functor2(mesh1, mesh2);
  const double distb1 = functor2(mesh2, mesh1);
  std::cout << std::endl;
  std::cout << "* Hausdorff distance1 f: " << dista0 << std::endl;
  std::cout << "* Hausdorff distance1 b: " << dista1 << std::endl;
  std::cout << std::endl;
  std::cout << "* Hausdorff distance2 f: " << distb0 << std::endl;
  std::cout << "* Hausdorff distance2 b: " << distb1 << std::endl;
}

template<
typename FunctionWrapper>
void test_timings(const std::string filepath, const FunctionWrapper& functor) {

  std::cout.precision(20);
  std::cout << std::endl << "* testing timing: " << functor.name() << std::endl;

  Timer timer;
  Surface_mesh mesh1, mesh2;
  get_meshes(filepath, filepath, mesh1, mesh2);

  timer.reset();
  timer.start();
  const double dista = functor(mesh1, mesh2);
  timer.stop();
  std::cout << "* time 0 (sec.): " << timer.time() << std::endl;

  PMP::transform(Affine_transformation_3(CGAL::Translation(),
    Vector_3(FT(0), FT(0), FT(1))), mesh2);

  timer.reset();
  timer.start();
  const double distb = functor(mesh1, mesh2);
  timer.stop();
  std::cout << "* time 1 (sec.): " << timer.time() << std::endl;

  std::cout << "dista = " << dista << std::endl;
  std::cout << "distb = " << distb << std::endl;
}

template<
typename FunctionWrapper>
void test_bunny(
  const FunctionWrapper& functor, const int n = 5, const bool save = false) {

  std::cout.precision(20);
  const std::string filepath1 = "data/bunny1.off"; // approx 16.3K
  const std::string filepath2 = "data/bunny2.off"; // approx 69.4K

  std::cout << std::endl << "* testing bunny:" << std::endl;
  std::cout << "* name -> " << functor.name() << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath1, filepath2, mesh1, mesh2);
  if (save) save_mesh(mesh1, "mesh1");
  if (save) save_mesh(mesh2, "mesh2");

  // Get 3 times longest dimension.
  const auto bbox = PMP::bbox(mesh2);
  const FT dist1 = static_cast<FT>(CGAL::abs(bbox.xmax() - bbox.xmin()));
  const FT dist2 = static_cast<FT>(CGAL::abs(bbox.ymax() - bbox.ymin()));
  const FT dist3 = static_cast<FT>(CGAL::abs(bbox.zmax() - bbox.zmin()));
  const FT dim = FT(3) * (CGAL::max)((CGAL::max)(dist1, dist2), dist3);

  // Get timings.
  Timer timer;
  std::vector<double> times;
  times.reserve(n);

  if (n == 0) {

    const FT distance = dim;
    PMP::transform(Affine_transformation_3(CGAL::Translation(),
      Vector_3(distance, distance, distance)), mesh2);
    if (save) save_mesh(mesh2, "mesh2");
    timer.reset();
    timer.start();
    const double dista = functor(mesh1, mesh2);
    const double distb = functor(mesh2, mesh1);
    const double distc = (CGAL::max)(dista, distb);
    timer.stop();
    times.push_back(timer.time());
    std::cout << "* distance / Hausdorff / time (sec.) : " <<
      distance << " / " << distc << " / " << times.back() << std::endl;

  } else {

    // t is the step where n is the number of steps.
    const FT t = FT(1) / static_cast<FT>(n);
    for (int k = n; k >= 0; --k) {
      auto mesh = mesh2;
      const FT distance = k * t * dim;
      PMP::transform(Affine_transformation_3(CGAL::Translation(),
        Vector_3(distance, distance, distance)), mesh);
      if (save) save_mesh(mesh, "mesh2-" + std::to_string(k));

      timer.reset();
      timer.start();
      const double dista = functor(mesh1, mesh);
      const double distb = functor(mesh, mesh1);
      const double distc = (CGAL::max)(dista, distb);
      timer.stop();
      times.push_back(timer.time());
      std::cout << "* distance / Hausdorff / time (sec.) " << k << " : " <<
        distance << " / " << distc << " / " << times.back() << std::endl;
    }
  }

  double min_time = +std::numeric_limits<double>::infinity();
  double max_time = -std::numeric_limits<double>::infinity();
  double avg_time = 0.0;
  for (const double time : times) {
    min_time = (CGAL::min)(min_time, time);
    max_time = (CGAL::max)(max_time, time);
    avg_time += time;
  }
  avg_time /= static_cast<double>(times.size());

  std::cout << std::endl << "* timings (msec.): " << std::endl;
  std::cout << "avg time: " << avg_time * 1000.0 << std::endl;
  std::cout << "min time: " << min_time * 1000.0 << std::endl;
  std::cout << "max time: " << max_time * 1000.0 << std::endl;
}

int main(int argc, char** argv) {

  // std::string name;
  // std::cin >> name;

  const double error_bound = 1e-4;
  const std::size_t num_samples = 1000;
  std::cout << std::endl << "* error bound: " << error_bound << std::endl;
  std::cout << std::endl << "* number of samples: " << num_samples << std::endl;
  std::string filepath = (argc > 1 ? argv[1] : "data/blobby.off");

  // ------------------------------------------------------------------------ //
  // Examples.

  // remeshing_tetrahedon_example(error_bound);
  // interior_triangle_example(error_bound);
  // perturbing_mesh_example(filepath, error_bound);
  // moving_mesh_example(filepath, filepath, 5, error_bound);

  // ------------------------------------------------------------------------ //
  // Tests.

  Approximate_hd_wrapper apprx_hd(num_samples);
  Naive_bounded_error_hd_wrapper naive_hd(error_bound);
  Bounded_error_hd_wrapper bound_hd(error_bound);

  // --- Testing basic properties.

  // test_synthetic_data(apprx_hd);
  // test_synthetic_data(naive_hd);
  test_synthetic_data(bound_hd);

  // --- Compare on common meshes.

  // test_one_versus_another(apprx_hd, naive_hd);
  // test_one_versus_another(naive_hd, bound_hd);
  test_one_versus_another(bound_hd, apprx_hd);

  // --- Compare on real meshes.

  const std::string filepath1 = (argc > 1 ? argv[1] : "data/blobby.off");
  const std::string filepath2 = (argc > 2 ? argv[2] : "data/tetrahedron-remeshed.off");

  // test_real_meshes(filepath1, filepath2, apprx_hd, naive_hd);
  // test_real_meshes(filepath1, filepath2, naive_hd, bound_hd);
  test_real_meshes(filepath1, filepath2, bound_hd, apprx_hd);

  // --- Compare timings.

  filepath = (argc > 1 ? argv[1] : "data/blobby-remeshed.off");
  // test_timings(filepath, apprx_hd);
  // test_timings(filepath, naive_hd);
  test_timings(filepath, bound_hd);

  // --- Compare with the paper.

  // test_bunny(apprx_hd);
  // test_bunny(naive_hd);
  test_bunny(bound_hd, 0);

  // ------------------------------------------------------------------------ //
  std::cout << std::endl;
}
