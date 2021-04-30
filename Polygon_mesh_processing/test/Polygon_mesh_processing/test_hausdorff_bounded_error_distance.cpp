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

using Kernel     = EPICK;
using FT         = typename Kernel::FT;
using Point_3    = typename Kernel::Point_3;
using Vector_3   = typename Kernel::Vector_3;
using Triangle_3 = typename Kernel::Triangle_3;

using TAG                     = CGAL::Parallel_if_available_tag;
using Surface_mesh            = CGAL::Surface_mesh<Point_3>;
using Affine_transformation_3 = CGAL::Aff_transformation_3<Kernel>;
using Timer                   = CGAL::Real_timer;

using Face_handle = typename boost::graph_traits<Surface_mesh>::face_descriptor;
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
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

template<class FunctionWrapper>
void test_synthetic_data(const FunctionWrapper& functor) {

  std::cout << std::endl << "-- test synthetic data:" << std::endl << std::endl;
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

  std::cout << std::endl << "-- test one versus another (tetrahedron):" << std::endl << std::endl;
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
  std::cout << std::endl << "-- test real meshes:" << std::endl << std::endl;
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

template<class FunctionWrapper>
void test_timings(const std::string filepath, const FunctionWrapper& functor) {

  std::cout.precision(20);
  std::cout << std::endl << "-- test timings: " << functor.name() << std::endl << std::endl;

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

  std::cout << "* dista = " << dista << std::endl;
  std::cout << "* distb = " << distb << std::endl;
}

template<class FunctionWrapper>
void test_bunny(
  const FunctionWrapper& functor, const int n = 5, const bool save = false) {

  std::cout.precision(20);
  const std::string filepath1 = "data/bunny1.off"; // approx 16.3K
  const std::string filepath2 = "data/bunny2.off"; // approx 69.4K

  std::cout << std::endl << "-- test bunny:" << std::endl << std::endl;
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
    // const double dista = functor(mesh1, mesh2);
    const double distb = functor(mesh2, mesh1);
    const double distc = distb; // (CGAL::max)(dista, distb);
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
  std::cout << "* avg time: " << avg_time * 1000.0 << std::endl;
  std::cout << "* min time: " << min_time * 1000.0 << std::endl;
  std::cout << "* max time: " << max_time * 1000.0 << std::endl;
}

Triangle_3 get_triangle(const int index, const Surface_mesh& mesh) {

  const auto mfaces = faces(mesh);
  assert(index >= 0);
  assert(index < static_cast<int>(mfaces.size()));
  const auto face = *(mfaces.begin() + index);

  const auto he = halfedge(face, mesh);
  const auto vertices = vertices_around_face(he, mesh);
  auto vit = vertices.begin();
  const auto& p0 = mesh.point(*vit); ++vit;
  const auto& p1 = mesh.point(*vit); ++vit;
  const auto& p2 = mesh.point(*vit);
  return Triangle_3(p0, p1, p2);
}

void compute_realizing_triangles(
  const Surface_mesh& mesh1, const Surface_mesh& mesh2,
  const double error_bound, const std::string prefix, const bool save = false) {

  std::cout << "* getting realizing triangles: " << std::endl;
  const double hdist =
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound);
  const auto data = std::make_pair(hdist, std::make_pair(Face_handle(), Face_handle()));
  const auto& faces = data.second;
  const int f1 = static_cast<int>(faces.first);
  const int f2 = static_cast<int>(faces.second);

  std::cout << "* Hausdorff: " << data.first << std::endl;
  std::cout << "* f1 / f2: " << f1 << " / " << f2 << std::endl;

  if (f1 != -1 && save) {
    const auto triangle = get_triangle(f1, mesh1);
    Surface_mesh sm1;
    sm1.add_vertex(triangle[0]);
    sm1.add_vertex(triangle[1]);
    sm1.add_vertex(triangle[2]);
    sm1.add_face(sm1.vertices());
    save_mesh(sm1, prefix + "triangle1");
  }

  if (f2 != -1 && save) {
    const auto triangle = get_triangle(f2, mesh2);
    Surface_mesh sm2;
    sm2.add_vertex(triangle[0]);
    sm2.add_vertex(triangle[1]);
    sm2.add_vertex(triangle[2]);
    sm2.add_face(sm2.vertices());
    save_mesh(sm2, prefix + "triangle2");
  }
}

void test_realizing_triangles(
  const double error_bound, const bool save = false) {

  std::cout.precision(20);
  std::cout << std::endl << "-- test realizing triangles:" << std::endl << std::endl;

  // Basic test.
  std::cout << " ---- basic test ---- " << std::endl;
  Surface_mesh mesh1, mesh2;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  mesh2.add_vertex(Point_3(0, 0, error_bound / 2.0));
  mesh2.add_vertex(Point_3(2, 0, error_bound / 2.0));
  mesh2.add_vertex(Point_3(1, 1, error_bound / 2.0));
  mesh2.add_face(mesh2.vertices());

  if (save) save_mesh(mesh1, "1mesh1");
  if (save) save_mesh(mesh2, "1mesh2");

  compute_realizing_triangles(mesh1, mesh2, error_bound, "1", save);

  mesh2.clear();
  mesh2.add_vertex(Point_3(0, 0, 1));
  mesh2.add_vertex(Point_3(2, 0, 1));
  mesh2.add_vertex(Point_3(1, 1, 1));
  mesh2.add_face(mesh2.vertices());

  if (save) save_mesh(mesh2, "1mesh2");
  compute_realizing_triangles(mesh1, mesh2, error_bound, "1", save);

  // Complex test.
  std::cout << std::endl << " ---- complex test ---- " << std::endl;
  const std::string filepath1 = "data/tetrahedron.off";
  const std::string filepath2 = "data/tetrahedron-remeshed.off";
  get_meshes(filepath1, filepath2, mesh1, mesh2);

  PMP::transform(Affine_transformation_3(
    CGAL::Translation(), Vector_3(0, 0, error_bound / 2.0)), mesh2);

  if (save) save_mesh(mesh1, "2mesh1");
  if (save) save_mesh(mesh2, "2mesh2");

  compute_realizing_triangles(mesh1, mesh2, error_bound, "2", save);

  mesh2.clear();
  get_mesh(filepath2, mesh2);
  PMP::transform(Affine_transformation_3(CGAL::Translation(),
    Vector_3(FT(0), FT(0), FT(1))), mesh2);

  if (save) save_mesh(mesh2, "2mesh2");
  compute_realizing_triangles(mesh1, mesh2, error_bound, "2", save);
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
  // Tests.

  // Approximate_hd_wrapper does not work with EPECK!
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
  // filepath = "/Users/monet/Documents/fork/pull-requests/hausdorff/data/bunny-dense.off";
  // test_timings(filepath, apprx_hd);
  // test_timings(filepath, naive_hd);
  test_timings(filepath, bound_hd);


  // --- Compare with the paper.

  // test_bunny(apprx_hd);
  // test_bunny(naive_hd);
  test_bunny(bound_hd, 3);


  // --- Test realizing triangles.
  // test_realizing_triangles(error_bound);


  // ------------------------------------------------------------------------ //

  std::cout << std::endl;
  return EXIT_SUCCESS;
}
