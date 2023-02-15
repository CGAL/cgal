// Use it to include parallel computations in the bounded error Hausdorff distance.
#define USE_PARALLEL_BEHD

// Use this def in order to get all DEBUG info related to the bounded-error Hausdorff code!
#define CGAL_HAUSDORFF_DEBUG

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Real_timer.h>

using SCK   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel     = EPICK;
using FT         = typename Kernel::FT;
using Point_3    = typename Kernel::Point_3;
using Vector_3   = typename Kernel::Vector_3;
using Triangle_3 = typename Kernel::Triangle_3;

using TAG                     = CGAL::Sequential_tag;
using Surface_mesh            = CGAL::Surface_mesh<Point_3>;
using Polyhedron              = CGAL::Surface_mesh<Point_3>;
using Affine_transformation_3 = CGAL::Aff_transformation_3<Kernel>;
using Timer                   = CGAL::Real_timer;

using Face_handle = typename boost::graph_traits<Surface_mesh>::face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

struct Approximate_hd_wrapper
{
  const std::size_t m_num_samples;
  const double m_error_bound;

  std::string name() const { return "approximate"; }

  explicit Approximate_hd_wrapper(const std::size_t num_samples,
                                  const double  error_bound)
    : m_num_samples(num_samples), m_error_bound(error_bound)
  { }

  double operator()(const Surface_mesh& tm1, const Surface_mesh& tm2) const
  {
    return PMP::approximate_Hausdorff_distance<TAG>(tm1, tm2,
      CGAL::parameters::number_of_points_on_faces(m_num_samples).random_seed(0),
      CGAL::parameters::number_of_points_on_faces(m_num_samples).random_seed(0));
  }

  double symmetric(const Surface_mesh& tm1, const Surface_mesh& tm2) const
  {
    return PMP::approximate_symmetric_Hausdorff_distance<TAG>(tm1, tm2,
      CGAL::parameters::number_of_points_on_faces(m_num_samples).random_seed(0),
      CGAL::parameters::number_of_points_on_faces(m_num_samples).random_seed(0));
  }

  bool check_value(const double val, const double expected) const
  {
    const bool res = ((CGAL::abs(val - expected)) <= m_error_bound);
    if(!res)
      std::cerr << "Expected " << expected << " but got " << val << std::endl;

    return res;
  }
};

struct Naive_bounded_error_hd_wrapper
{
  const double m_error_bound;

  std::string name() const { return "naive bounded error"; }

  explicit Naive_bounded_error_hd_wrapper(const double error_bound) : m_error_bound(error_bound) { }

  double operator()(const Surface_mesh& tm1, const Surface_mesh& tm2) const
  {
    return PMP::bounded_error_Hausdorff_distance_naive<TAG>(tm1, tm2, m_error_bound);
  }

  double symmetric(const Surface_mesh& tm1, const Surface_mesh& tm2) const
  {
    const double dista = operator()(tm1, tm2);
    const double distb = operator()(tm2, tm1);
    return (CGAL::max)(dista, distb);
  }

  bool check_value(const double val, const double expected) const
  {
    const bool res = ((CGAL::abs(val - expected)) <= m_error_bound);
    if(!res)
      std::cerr << "Expected " << expected << " but got " << val << std::endl;

    return res;
  }
};

struct Bounded_error_hd_wrapper
{
  const double m_error_bound;

  std::string name() const { return "bounded error"; }

  explicit Bounded_error_hd_wrapper(const double error_bound) : m_error_bound(error_bound) { }

  double operator()(const Surface_mesh& tm1, const Surface_mesh& tm2) const
  {
    return PMP::bounded_error_Hausdorff_distance<TAG>(tm1, tm2, m_error_bound);
  }

  double symmetric(const Surface_mesh& tm1, const Surface_mesh& tm2) const
  {
    return PMP::bounded_error_symmetric_Hausdorff_distance<TAG>(tm1, tm2, m_error_bound);
  }

  bool check_value(const double val, const double expected) const
  {
    const bool res = ((CGAL::abs(val - expected)) <= m_error_bound);
    if(!res)
      std::cerr << "Expected " << expected << " but got " << val << std::endl;

    return res;
  }
};

template<typename PolygonMesh>
void get_mesh(const std::string& filepath, PolygonMesh& mesh)
{
  mesh.clear();
  CGAL::IO::read_polygon_mesh(filepath, mesh);
  std::cout << "* getting mesh with " << faces(mesh).size() << " faces" << std::endl;
}

template<typename PolygonMesh1, typename PolygonMesh2>
void get_meshes(const std::string& filepath1,
                const std::string& filepath2,
                PolygonMesh1& mesh1,
                PolygonMesh2& mesh2)
{
  get_mesh(filepath1, mesh1);
  get_mesh(filepath2, mesh2);
}

template<typename PolygonMesh>
void save_mesh(const PolygonMesh& mesh,
               const std::string& filepath)
{
  if(!CGAL::IO::write_PLY(filepath + ".ply", mesh))
  {
    std::cerr << "ERROR: cannot save this file: " << filepath << std::endl;
    exit(EXIT_FAILURE);
  }
}

// An easy example of a tetrahedron and its remeshed version.
void remeshing_tetrahedon_example(const double error_bound,
                                  const bool save = true)
{
  Timer timer;
  Surface_mesh mesh1, mesh2;
  std::cout << std::endl << "* (E1) remeshing Tetrahedron example:" << std::endl;

  CGAL::make_tetrahedron(Point_3(0, 0, 0), Point_3(2, 0, 0),
                         Point_3(1, 1, 1), Point_3(1, 0, 2), mesh1);
  mesh2 = mesh1;

  using edge_descriptor = typename boost::graph_traits<Surface_mesh>::edge_descriptor;
  Surface_mesh::Property_map<edge_descriptor, bool> is_constrained_map =
    mesh2.add_property_map<edge_descriptor, bool>("e:is_constrained", true).first;

  const double target_edge_length = 0.1;
  PMP::isotropic_remeshing(mesh2.faces(), target_edge_length, mesh2,
                           CGAL::parameters::edge_is_constrained_map(is_constrained_map));

  std::cout << "Mesh1: " << num_vertices(mesh1) << " nv " << num_faces(mesh1) << " nf" << std::endl;
  std::cout << "Mesh2: " << vertices(mesh2).size() << " nv " << faces(mesh2).size() << " nf" << std::endl;

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  timer.reset();
  timer.start();
  const double hdist = PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound);
  std::cout << "* bounded-error Hausdorff distance: " << hdist << std::endl;
  timer.stop();
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;
  assert(hdist <= error_bound); // real theoretical distance is 0
}

// Example with a point realizing the Hausdorff distance strictly lying in the interior of a triangle.
void interior_triangle_example(const double error_bound,
                               const bool save = true)
{
  Timer timer;
  Surface_mesh mesh1, mesh2;
  std::cout << std::endl << "* (E2) interior Triangle example:" << std::endl;

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
  mesh2.add_face(v3, v4, v5);
  mesh2.add_face(v0, v4, v3);
  mesh2.add_face(v1, v5, v4);
  mesh2.add_face(v2, v3, v5);

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  timer.reset();
  timer.start();
  const double hdist = PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound);
  timer.stop();
  std::cout << "* bounded-error Hausdorff distance: " << hdist << std::endl;
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;

  assert(1.59168 < hdist && hdist < 1.59169);
}

// Read a real mesh given by the user, perturb it slightly, and compute the
// Hausdorff distance between the original mesh and its perturbation.
void perturbing_surface_mesh_example(const std::string& filepath,
                                     const double error_bound,
                                     const bool save = true)
{
  std::cout << std::endl << "* (E3) perturbing Surface Mesh example:" << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath, filepath, mesh1, mesh2);

  const double max_size = 0.1;
  std::cout << "* perturbing the second mesh; max_size: " << max_size << std::endl;
  PMP::random_perturbation(mesh2.vertices(), mesh2, max_size, CGAL::parameters::do_project(false)
                                                                               .random_seed(0));

  if(save)
  {
    save_mesh(mesh2, "mesh1");
    save_mesh(mesh1, "mesh2");
  }

  Timer timer;
  timer.start();
  const double hdist = PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound);
  timer.stop();
  std::cout << "* bounded-error Hausdorff distance: " << hdist << std::endl;
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;

  assert(hdist > 0. && hdist <= max_size);
}

// Read two meshes and store them in two different face graph containers,
// perturb the second mesh, and compute the Hausdorff distance.
void perturbing_polyhedron_mesh_example(const std::string& filepath,
                                        const double error_bound,
                                        const bool save = true)
{
  std::cout << std::endl << "* (E3 bis) perturbing Polyhedron mesh example:" << std::endl;

  Surface_mesh mesh1;
  Polyhedron mesh2;
  get_mesh(filepath, mesh1);
  get_mesh(filepath, mesh2);

  const double max_size = 0.1;
  std::cout << "* perturbing the second mesh; max_size: " << max_size << std::endl;
  PMP::random_perturbation(vertices(mesh2), mesh2, max_size, CGAL::parameters::do_project(false)
                                                                              .random_seed(0));

  if(save)
  {
    save_mesh(mesh2, "mesh1");
    save_mesh(mesh1, "mesh2");
  }

  Timer timer;
  timer.start();
  const double hdist1 = PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound);
  std::cout << "* bounded-error Hausdorff distance 1->2: " << hdist1 << std::endl;
  const double hdist2 = PMP::bounded_error_Hausdorff_distance<TAG>(mesh2, mesh1, error_bound);
  std::cout << "* bounded-error Hausdorff distance 2->1: " << hdist2 << std::endl;
  timer.stop();
  std::cout << "* processing time: " << timer.time() << " s." << std::endl;
  assert(hdist1 > 0. && hdist1 <= max_size);
  assert(hdist2 > 0.);
  assert(hdist2 > hdist1);

  const double hdist3 = PMP::bounded_error_symmetric_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound);
  std::cout << "* symmetric bounded-error Hausdorff distance: " << hdist3 << std::endl;
  assert(hdist3 == (CGAL::max)(hdist1, hdist2));
}

// Read two meshes given by the user, initially place them at their originally
// given position. Move the second mesh in 300 steps away from the first one.
// Print how the Hausdorff distance changes.
void moving_surface_mesh_example(const std::string& filepath1,
                                 const std::string& filepath2,
                                 const std::size_t n,
                                 const double error_bound,
                                 const bool save = true)
{
  Timer timer;
  std::cout << std::endl << "* (E4) moving Surface Mesh example:" << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath1, filepath2, mesh1, mesh2);

  const CGAL::Bbox_3 bbox = PMP::bbox(mesh2);
  const FT distance = CGAL::approximate_sqrt(CGAL::squared_distance(Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin()),
                                                                    Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax())));

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const FT t = FT(1) / FT(100);
  const FT a = t * distance;
  const Vector_3 v(a, a, a);
  const FT vl = std::sqrt(3) * t * distance;

  for(std::size_t i=0; i<n; ++i)
  {
    PMP::transform(Affine_transformation_3(CGAL::Translation(), v), mesh2);
    if(save)
      save_mesh(mesh2, "mesh2-" + std::to_string(i + 1));

    timer.reset();
    timer.start();
    const double hdist = PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound);
    timer.stop();
    std::cout << "* position: " << i << " translation distance: " << (i+1) * vl << std::endl;
    std::cout << "* bounded-error Hausdorff distance: " << hdist << std::endl;
    std::cout << "* processing time: " << timer.time() << " s." << std::endl;
    const FT l = (i+1) * vl - error_bound, u = (i+1) * vl + error_bound;
    assert(hdist >= l && hdist <= u);
  }
}

template<class FunctionWrapper>
void test_0(const FunctionWrapper& functor,
            const bool save = true)
{
  // The same triangle.
  // Expected distance is 0.

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 0 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());
  mesh2 = mesh1;

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 0): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 0): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 0): " << distc << std::endl;

  assert(functor.check_value(dista, 0));
  assert(functor.check_value(distb, 0));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_1(const FunctionWrapper& functor,
            const bool save = true)
{
  // Two parallel triangles, 1 unit distance away from each other.
  // Expected distance is 1.

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 1 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  mesh2.add_vertex(Point_3(0, 0, 1));
  mesh2.add_vertex(Point_3(2, 0, 1));
  mesh2.add_vertex(Point_3(1, 1, 1));
  mesh2.add_face(mesh2.vertices());

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 1): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 1): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 1): " << distc << std::endl;

  assert(functor.check_value(dista, 1.));
  assert(functor.check_value(distb, 1.));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_2(const FunctionWrapper& functor,
            const bool save = true)
{
  // One triangle is orthogonal to the other one and shares a common edge.
  // Expected distance is 1.

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 2 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  mesh2.add_vertex(Point_3(0, 0, 0));
  mesh2.add_vertex(Point_3(2, 0, 0));
  mesh2.add_vertex(Point_3(1, 0, 1));
  mesh2.add_face(mesh2.vertices());

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 1): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 1): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 1): " << distc << std::endl;

  assert(functor.check_value(dista, 1.));
  assert(functor.check_value(distb, 1.));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_3(const FunctionWrapper& functor,
            const bool save = true)
{
  // One triangle is orthogonal to the other one and shares a common edge
  // that is moved 1 unit distance away.
  // Expected distances are 2 and sqrt(2) ~= 1.4142135623730951.

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 3 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  mesh2.add_vertex(Point_3(0, 0, 1));
  mesh2.add_vertex(Point_3(2, 0, 1));
  mesh2.add_vertex(Point_3(1, 0, 2));
  mesh2.add_face(mesh2.vertices());

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 1.4142135623730951): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 2): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 2): " << distc << std::endl;

  assert(functor.check_value(dista, std::sqrt(2.)));
  assert(functor.check_value(distb, 2.));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_4(const FunctionWrapper& functor,
            const bool save = true)
{
  // One triangle is orthogonal to the other one and shares a common vertex.
  // Expected distance is std::sqrt(2) * std::sin(60°) ~= 1.2247448713915892

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 4 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  mesh2.add_vertex(Point_3(0, 1, 1));
  mesh2.add_vertex(Point_3(2, 1, 1));
  mesh2.add_vertex(Point_3(1, 1, 0));
  mesh2.add_face(mesh2.vertices());

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 1.2247448713915892): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 1.2247448713915892): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 1.2247448713915892): " << distc << std::endl;

  assert(functor.check_value(dista, 1.2247448713915892));
  assert(functor.check_value(distb, 1.2247448713915892));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_5(const FunctionWrapper& functor, const bool save = true) {

  // One triangle is orthogonal to the other one and shares a common vertex
  // that is moved 1 unit distance away.
  // Expected distances are sqrt(3) ~= 1.7320508075688771 and sqrt(4.5) ~= 2.1213203435596423.

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 5 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  mesh2.add_vertex(Point_3(0, 1, 2));
  mesh2.add_vertex(Point_3(2, 1, 2));
  mesh2.add_vertex(Point_3(1, 1, 1));
  mesh2.add_face(mesh2.vertices());

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 1.7320508075688771): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 2.1213203435596423): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 2.1213203435596423): " << distc << std::endl;

  assert(functor.check_value(dista, 1.7320508075688771));
  assert(functor.check_value(distb, 2.1213203435596423));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_6(const FunctionWrapper& functor,
            const bool save = true)
{
  // The first and second mesh have different number of triangles.
  // They are parallel and lie in the same plane. The middle triangle is common to both meshes.
  // Expected distances are 0 and 1/sqrt(2) ~=0.70710678118654746.

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 6 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  const auto v0 = mesh2.add_vertex(Point_3(0, 0, 0));
  const auto v1 = mesh2.add_vertex(Point_3(2, 0, 0));
  const auto v2 = mesh2.add_vertex(Point_3(1, 1, 0));
  const auto v3 = mesh2.add_vertex(Point_3(2, 1, 0));
  const auto v4 = mesh2.add_vertex(Point_3(0, 1, 0));

  mesh2.add_face(v0, v1, v2);
  mesh2.add_face(v2, v1, v3);
  mesh2.add_face(v0, v2, v4);

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 0): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 0.70710678118654746): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 0.70710678118654746): " << distc << std::endl;

  assert(functor.check_value(dista, 0.));
  assert(functor.check_value(distb, 0.70710678118654746));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_7(const FunctionWrapper& functor,
            const bool save = true)
{
  // One triangle is moved to 0.5 unit distance away from 3 other triangles.
  // The first and second meshes are parallel.
  // Expected distances are 0.5 and sqrt(0.75) ~= 0.8660254037844386.

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 7 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  const auto v0 = mesh2.add_vertex(Point_3(0, 0, 0.5));
  const auto v1 = mesh2.add_vertex(Point_3(2, 0, 0.5));
  const auto v2 = mesh2.add_vertex(Point_3(1, 1, 0.5));
  const auto v3 = mesh2.add_vertex(Point_3(2, 1, 0.5));
  const auto v4 = mesh2.add_vertex(Point_3(0, 1, 0.5));
  mesh2.add_face(v0, v1, v2);
  mesh2.add_face(v2, v1, v3);
  mesh2.add_face(v0, v2, v4);

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 0.5): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 0.8660254037844386): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 0.8660254037844386): " << distc << std::endl;

  assert(functor.check_value(dista, 0.5));
  assert(functor.check_value(distb, 0.8660254037844386));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_8(const FunctionWrapper& functor,
            const bool save = true)
{
  // One mesh has one triangle at zero level, another mesh has two triangles
  // where the first one is at level 1 and the second one is at level 2.
  // Expected distances are 1 and 2.

  Surface_mesh mesh1, mesh2;
  std::cout << " ---- testing 8 ---- " << std::endl;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  const auto v0 = mesh2.add_vertex(Point_3(0, 0, 1));
  const auto v1 = mesh2.add_vertex(Point_3(2, 0, 1));
  const auto v2 = mesh2.add_vertex(Point_3(1, 1, 1));
  const auto v3 = mesh2.add_vertex(Point_3(0, 0, 2));
  const auto v4 = mesh2.add_vertex(Point_3(2, 0, 2));
  const auto v5 = mesh2.add_vertex(Point_3(1, 1, 2));
  mesh2.add_face(v0, v1, v2);
  mesh2.add_face(v3, v4, v5);

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 1): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 2): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 2): " << distc << std::endl;

  assert(functor.check_value(dista, 1.));
  assert(functor.check_value(distb, 2.));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_9(const FunctionWrapper& functor,
            const bool save = true)
{
  // Two meshes partially overlap, have 2 triangles in common and each one has
  // two its own trianles. All triangles form a Z shape where the height is 1.
  // The expected result is 1.

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

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  const double dista = functor(mesh1, mesh2);
  const double distb = functor(mesh2, mesh1);
  const double disto = (CGAL::max)(dista, distb);
  const double distc = functor.symmetric(mesh1, mesh2);

  std::cout << "* Hausdorff distance (expected 1): " << dista << std::endl;
  std::cout << "* HInverted distance (expected 1): " << distb << std::endl;
  std::cout << "* Symmetric distance (expected 1): " << distc << std::endl;

  assert(functor.check_value(dista, 1.));
  assert(functor.check_value(distb, 1.));
  assert(functor.check_value(disto, distc));
}

template<class FunctionWrapper>
void test_synthetic_data(const FunctionWrapper& functor)
{
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

template <typename FunctionWrapper1,
          typename FunctionWrapper2>
void test_one_versus_another(const FunctionWrapper1& functor1,
                             const FunctionWrapper2& functor2,
                             const double error_bound)
{
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
  std::cout << "* Hausdorff distance1: " << dista0 << std::endl;
  const double distb0 = functor2(mesh1, mesh2);
  std::cout << "* Hausdorff distance2: " << distb0 << std::endl;
  assert(dista0 <= error_bound);
  assert(distb0 <= error_bound);

  std::cout << "* translating by 1 unit ..." << std::endl;
  PMP::transform(Affine_transformation_3(CGAL::Translation(), Vector_3(FT(0), FT(0), FT(1))), mesh2);

  // TEST 1.
  // Translate by 1 unit distance and compare.
  // The expected distance is 1.
  std::cout << " ---- testing 1 ---- " << std::endl;
  const double dista1 = functor1(mesh1, mesh2);
  std::cout << "* Hausdorff distance1: " << dista1 << std::endl;
  const double distb1 = functor2(mesh1, mesh2);
  std::cout << "* Hausdorff distance2: " << distb1 << std::endl;

  assert(std::abs(dista0 - distb0) <= 2. * error_bound);
  assert(std::abs(dista1 - distb1) <= 2. * error_bound);
}

template <typename FunctionWrapper1,
          typename FunctionWrapper2>
void test_real_meshes(const std::string& filepath1,
                      const std::string& filepath2,
                      const FunctionWrapper1& functor1,
                      const FunctionWrapper2& functor2,
                      const double error_bound)
{
  std::cout << std::endl << "-- test real meshes:" << std::endl << std::endl;
  std::cout << "* name 1 -> " << functor1.name() << std::endl;
  std::cout << "* name 2 -> " << functor2.name() << std::endl;
  std::cout << "* input path 1: " << filepath1 << std::endl;
  std::cout << "* input path 2: " << filepath2 << std::endl;

  Surface_mesh mesh1, mesh2;
  get_meshes(filepath1, filepath2, mesh1, mesh2);

  // Load and compare.
  auto run = [&](auto f, const Surface_mesh& m1, const Surface_mesh& m2) -> std::pair<double, double>
  {
    Timer timer;
    timer.start();
    const double dist = f(m1, m2);
    timer.stop();
    return std::make_pair(dist, timer.time());
  };

  double timea0, timeb0, timea1, timeb1;
  double dista0, distb0, dista1, distb1;
  std::tie(dista0, timea0) = run(functor1, mesh1, mesh2);
  std::tie(distb0, timeb0) = run(functor2, mesh1, mesh2);
  std::tie(dista1, timea1) = run(functor1, mesh2, mesh1);
  std::tie(distb1, timeb1) = run(functor2, mesh2, mesh1);

  std::cout << std::endl;
  std::cout << "* Hausdorff distance1-2 " << functor1.name() << " dist " << dista0 << " runtime: " << timea0 << std::endl;
  std::cout << "* Hausdorff distance1-2   " << functor2.name() << " dist " << distb0 << " runtime: " << timeb0 << std::endl;
  std::cout << std::endl;
  std::cout << "* Hausdorff distance2-1 " << functor1.name() << " dist " << dista1 << " runtime: " << timea1 << std::endl;
  std::cout << "* Hausdorff distance2-1   " << functor2.name() << " dist " << distb1 << " runtime: " << timeb1 << std::endl;

  assert(std::abs(dista0 - distb0) <= 2. * error_bound);
  assert(std::abs(dista1 - distb1) <=  2. * error_bound);
}

Triangle_3 get_triangle(const int index,
                        const Surface_mesh& mesh)
{
  const auto mfaces = faces(mesh);
  assert(index >= 0);
  assert(index < static_cast<int>(mfaces.size()));
  const auto face = *(std::next(mfaces.begin(), index));

  const auto he = halfedge(face, mesh);
  const auto vertices = vertices_around_face(he, mesh);
  assert(vertices.size() >= 3);

  auto vit = vertices.begin();
  const Point_3& p0 = mesh.point(*vit); ++vit;
  const Point_3& p1 = mesh.point(*vit); ++vit;
  const Point_3& p2 = mesh.point(*vit);

  return Triangle_3(p0, p1, p2);
}

void compute_realizing_triangles(const Surface_mesh& mesh1,
                                 const Surface_mesh& mesh2,
                                 const double error_bound,
                                 const int expected_f1,
                                 const int expected_f2,
                                 const std::string& prefix,
                                 const bool save = true)
{
  std::cout << "* getting realizing triangles: " << std::endl;
  std::vector< std::pair<Face_handle, Face_handle> > fpairs;
  const double hdist = PMP::bounded_error_Hausdorff_distance<TAG>(
                         mesh1, mesh2, error_bound,
                         CGAL::parameters::output_iterator(std::back_inserter(fpairs)));
  assert(fpairs.size() == 2); // lower bound face pair + upper bound face pair

  const int f1 = static_cast<int>(fpairs.back().first);  // f1 / f2: upper bound
  const int f2 = static_cast<int>(fpairs.back().second);

  std::cout << "* Hausdorff: " << hdist << std::endl;
  std::cout << "* f1 / f2: " << f1 << " / " << f2 << std::endl;
  std::cout << "* expected f1 / f2: " << expected_f1 << " / " << expected_f2 << std::endl;

  if(save)
  {
    auto triangle = get_triangle(f1, mesh1);
    Surface_mesh sm1;
    sm1.add_vertex(triangle[0]);
    sm1.add_vertex(triangle[1]);
    sm1.add_vertex(triangle[2]);
    sm1.add_face(sm1.vertices());
    save_mesh(sm1, prefix + "triangle1");

    triangle = get_triangle(f2, mesh2);
    Surface_mesh sm2;
    sm2.add_vertex(triangle[0]);
    sm2.add_vertex(triangle[1]);
    sm2.add_vertex(triangle[2]);
    sm2.add_face(sm2.vertices());
    save_mesh(sm2, prefix + "triangle2");
  }

  assert(f1 == expected_f1);
  assert(f2 == expected_f2);
}

void test_realizing_triangles(const double error_bound,
                              const bool save = true)
{
  std::cout << std::endl << "-- test realizing triangles:" << std::endl << std::endl;

  // Basic test.
  std::cout << " ---- basic test ---- " << std::endl;
  Surface_mesh mesh1, mesh2;

  mesh1.add_vertex(Point_3(0, 0, 0));
  mesh1.add_vertex(Point_3(2, 0, 0));
  mesh1.add_vertex(Point_3(1, 1, 0));
  mesh1.add_face(mesh1.vertices());

  mesh2.add_vertex(Point_3(0, 0, error_bound / 2.));
  mesh2.add_vertex(Point_3(2, 0, error_bound / 2.));
  mesh2.add_vertex(Point_3(1, 1, error_bound / 2.));
  mesh2.add_face(mesh2.vertices());

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  compute_realizing_triangles(mesh1, mesh2, error_bound, 0, 0, "1", save);

  std::cout << std::endl << " ---- complex test ---- " << std::endl;

  Surface_mesh tmp1,tmp2;
  std::array<Surface_mesh::Vertex_index,3> vhs1 =
  {
    tmp1.add_vertex(Point_3(0, 1, 3)),
    tmp1.add_vertex(Point_3(1, 1, 3)),
    tmp1.add_vertex(Point_3(1, 0, 3))
  };
  tmp1.add_face(vhs1);

  std::array<Surface_mesh::Vertex_index,3> vhs2 =
  {
    tmp2.add_vertex(Point_3(0, 1, 3.5)),
    tmp2.add_vertex(Point_3(1, 1, 3.5)),
    tmp2.add_vertex(Point_3(1, 0, 3.5))
  };
  tmp2.add_face(vhs2);

  const std::string filepath1 = "data/tetrahedron.off";
  const std::string filepath2 = "data/tetrahedron-remeshed.off";
  mesh1.clear();
  mesh2.clear();
  get_meshes(filepath1, filepath2, mesh1, mesh2);
  PMP::transform(Affine_transformation_3(CGAL::Translation(), Vector_3(0, 0, 10 * error_bound)), tmp2);

  mesh1.join(tmp1);
  mesh2.join(tmp2);

  if(save)
  {
    save_mesh(mesh1, "mesh1-2");
    save_mesh(mesh2, "mesh2-2");
  }

  const int expected_f1 = static_cast<int>(faces(mesh1).size() - 1);
  const int expected_f2 = static_cast<int>(faces(mesh2).size() - 1);
  compute_realizing_triangles(mesh1, mesh2, error_bound, expected_f1, expected_f2, "2", save);
}

#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED) && defined(USE_PARALLEL_BEHD)
void test_parallel_version(const std::string& filepath,
                           const double error_bound)
{
  std::cout << std::endl << "-- test parallel version:" << std::endl << std::endl;

  Timer timer;
  Surface_mesh mesh1, mesh2;
  get_meshes(filepath, filepath, mesh1, mesh2);

  // One-sided distance.
  std::cout << " ---- one-sided distance ---- " << std::endl;

  PMP::transform(Affine_transformation_3(CGAL::Translation(),
                                         Vector_3(FT(0), FT(0), FT(1))), mesh2);

  std::cout << " ---- SEQUENTIAL ---- " << std::endl;
  timer.reset();
  timer.start();
  const double dista = PMP::bounded_error_Hausdorff_distance<CGAL::Sequential_tag>(
    mesh1, mesh2, error_bound,
    CGAL::parameters::match_faces(false),
    CGAL::parameters::match_faces(false));
  timer.stop();
  const double timea = timer.time();

  std::cout << " ---- PARALLEL ---- " << std::endl;
  timer.reset();
  timer.start();
  const double distb = PMP::bounded_error_Hausdorff_distance<CGAL::Parallel_tag>(
    mesh1, mesh2, error_bound,
    CGAL::parameters::match_faces(false),
    CGAL::parameters::match_faces(false));
  timer.stop();
  const double timeb = timer.time();

  std::cout << "* time a seq (sec.): " << timea << std::endl;
  std::cout << "* time b par (sec.): " << timeb << std::endl;

  std::cout << "* dista seq = " << dista << std::endl;
  std::cout << "* distb par = " << distb << std::endl;

  assert(dista == distb);
}
#endif // defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED)

void test_early_quit(const std::string& filepath,
                     const bool save = true)
{
  std::cout << std::endl << "-- test early quit:" << std::endl << std::endl;

  Timer timer;
  Surface_mesh mesh1, mesh2;
  get_meshes(filepath, filepath, mesh1, mesh2);

  PMP::transform(Affine_transformation_3(CGAL::Translation(),
                                         Vector_3(0.0, 0.0, 0.5)), mesh2);

  if(save)
  {
    save_mesh(mesh1, "mesh1");
    save_mesh(mesh2, "mesh2");
  }

  std::cout << " ---- distance 0.5 < 1.0 ---- " << std::endl;
  timer.reset();
  timer.start();
  assert(!PMP::is_Hausdorff_distance_larger<TAG>(mesh1, mesh2, 1.0, 0.0001));
  timer.stop();
  const double timeb = timer.time();

  std::cout << std::endl << " ---- distance 0.5 < 0.6 ---- " << std::endl;
  timer.reset();
  timer.start();
  assert(!PMP::is_Hausdorff_distance_larger<TAG>(mesh1, mesh2, 0.6, 0.0001));
  timer.stop();
  const double timec = timer.time();

  std::cout << std::endl << " ---- distance 0.5 > 0.4 ---- " << std::endl;
  timer.reset();
  timer.start();
  assert(PMP::is_Hausdorff_distance_larger<TAG>(mesh1, mesh2, 0.4, 0.0001));
  timer.stop();
  const double timed = timer.time();

  std::cout << "* timeb 0.5 < 1.0; time: " << timeb << std::endl;
  std::cout << "* timec 0.5 < 0.6; time: " << timec << std::endl;
  std::cout << "* timed 0.5 > 0.4; time: " << timed << std::endl;
}

void run_examples(const double error_bound,
                  const std::string& filepath)
{
  remeshing_tetrahedon_example(error_bound);
  interior_triangle_example(error_bound);
  perturbing_surface_mesh_example(filepath, error_bound);
  perturbing_polyhedron_mesh_example(filepath, error_bound);
  moving_surface_mesh_example(filepath, filepath, 5, error_bound);
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const double error_bound = 1e-4;
  const std::size_t num_samples = 1000;
  std::cout << std::endl << "* error bound: " << error_bound << std::endl;
  std::cout << std::endl << "* number of samples: " << num_samples << std::endl;

  const std::string filepath = (argc > 1 ? argv[1] : CGAL::data_file_path("meshes/blobby.off"));
  run_examples(error_bound, filepath);

  // ------------------------------------------------------------------------ //
  // Tests.

  Approximate_hd_wrapper apprx_hd(num_samples, error_bound);
  Naive_bounded_error_hd_wrapper naive_hd(error_bound);
  Bounded_error_hd_wrapper bound_hd(error_bound);

  // --- Testing basic properties.
  test_synthetic_data(apprx_hd);
// test_synthetic_data(naive_hd); // naive takes too long
  test_synthetic_data(bound_hd);

  // --- Compare on common meshes.
//  test_one_versus_another(apprx_hd, naive_hd, error_bound); // naive takes too long
//  test_one_versus_another(naive_hd, bound_hd, error_bound); // naive takes too long
  test_one_versus_another(bound_hd, apprx_hd, error_bound);

  // --- Compare on real meshes.
  const std::string filepath1 = (argc > 1 ? argv[1] : CGAL::data_file_path("meshes/blobby.off"));
  const std::string filepath2 = (argc > 2 ? argv[2] : "data/tetrahedron-remeshed.off");

//  test_real_meshes(filepath1, filepath2, apprx_hd, naive_hd, error_bound);
//  test_real_meshes(filepath1, filepath2, naive_hd, bound_hd, error_bound);
  test_real_meshes(filepath1, filepath2, bound_hd, apprx_hd, error_bound);
//  test_real_meshes("data/elephant_concave_hole.off", CGAL::data_file_path("meshes/mech-holes-shark.off"), bound_hd, apprx_hd, error_bound); // commenting because approx hausdorff is annoyingly rough
  test_real_meshes("data/small_spheres.off", "data/overlapping_triangles.off", bound_hd, apprx_hd, error_bound);

  // --- Test realizing triangles.
  test_realizing_triangles(error_bound);

#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED) && defined(USE_PARALLEL_BEHD)
  // --- Test parallelization.
  test_parallel_version(filepath, error_bound);
#endif // defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED)

  // --- Test early quit.
  test_early_quit(filepath);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
