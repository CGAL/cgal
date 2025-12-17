#define CGAL_USE_OPTI_WITH_BBOX

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/internal/kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_integer.h>

#include <CGAL/Real_timer.h>

using K=CGAL::Exact_predicates_inexact_constructions_kernel;
using EK=CGAL::Exact_predicates_exact_constructions_kernel;
using Mesh=CGAL::Surface_mesh<K::Point_3>;
using EMesh=CGAL::Surface_mesh<EK::Point_3>;

namespace PMP=CGAL::Polygon_mesh_processing;

template<class Mesh>
void round_mesh(Mesh &m){
  auto exp = [](const double v)
  {
    int n;
    frexp(v, &n);
    return n;
  };
  auto pow_2 = [](const int k)
  {
      return (k>=0)?std::pow(2,k):1./std::pow(2,-k);
  };

  CGAL::Bbox_3 bb;
  for(auto vh : m.vertices()){
    bb += m.point(vh).bbox();
  }
  size_t grid_size=26;
  std::array<double, 3> max_abs{(std::max)(-bb.xmin(), bb.xmax()),
                                (std::max)(-bb.ymin(), bb.ymax()),
                                (std::max)(-bb.zmin(), bb.zmax())};
  // Compute scale so that the exponent of max absolute value are 52-1.
  std::array<double, 3> scale{pow_2(grid_size - exp(max_abs[0]) - 1),
                              pow_2(grid_size - exp(max_abs[1]) - 1),
                              pow_2(grid_size - exp(max_abs[2]) - 1)};

  auto round=[&](typename Mesh::Point &p){
    return typename Mesh::Point(std::ceil((CGAL::to_double(p.x()) * scale[0]) - 0.5) / scale[0],
                                std::ceil((CGAL::to_double(p.y()) * scale[1]) - 0.5) / scale[1],
                                std::ceil((CGAL::to_double(p.z()) * scale[2]) - 0.5) / scale[2]);
  };

  for(auto vh : m.vertices()){
    auto &p = m.point(vh);
    p = round(p);
  }
}

template<class Mesh>
void to_integer_mesh(Mesh &m){
  auto exp = [](const double v)
  {
    int n;
    frexp(v, &n);
    return n;
  };
  auto pow_2 = [](const int k)
  {
      return (k>=0)?std::pow(2,k):1./std::pow(2,-k);
  };

  CGAL::Bbox_3 bb;
  for(auto vh : m.vertices()){
    bb += m.point(vh).bbox();
  }
  size_t grid_size=26;
  std::array<double, 3> max_abs{(std::max)(-bb.xmin(), bb.xmax()),
                                (std::max)(-bb.ymin(), bb.ymax()),
                                (std::max)(-bb.zmin(), bb.zmax())};
  // Compute scale so that the exponent of max absolute value are 52-1.
  std::array<double, 3> scale{pow_2(grid_size - exp(max_abs[0]) - 1),
                              pow_2(grid_size - exp(max_abs[1]) - 1),
                              pow_2(grid_size - exp(max_abs[2]) - 1)};

  // auto round=[&](Point_3 &p){
  auto round=[&](typename Mesh::Point &p){
    return typename Mesh::Point(std::ceil((CGAL::to_double(p.x()) * scale[0]) - 0.5),
                                std::ceil((CGAL::to_double(p.y()) * scale[1]) - 0.5),
                                std::ceil((CGAL::to_double(p.z()) * scale[2]) - 0.5));
  };

  for(auto vh : m.vertices()){
    auto &p = m.point(vh);
    p = round(p);
  }
}

void test_traits()
{
  Mesh m;
  std::ifstream(CGAL::data_file_path("meshes/elephant.off")) >> m;

  std::pair p(0,0.05);
  using Traits = PMP::Orthogonal_cut_plane_traits<K>;

  PMP::clip(m, p, CGAL::parameters::geom_traits(Traits()));

  std::ofstream("clipped.off") << m;
}

void elementary_test_kernel()
{
  using Point_3 = typename EK::Point_3;
  CGAL::Real_timer timer;
  EMesh m;
  // make_hexahedron(Point_3(0, 1, 0), Point_3(0, 1, 1), Point_3(1, 0, 1), Point_3(1, 0, 0),
  //                 Point_3(0.75, 0.25, 0.25), Point_3(0.25, 0.75, 0.25), Point_3(0.25, 0.75, 0.75), Point_3(0.75, 0.25, 0.75),
  //                 m);
  // make_tetrahedron(Point_3(1, 0, 0), Point_3(0, 1, 0), Point_3(0, 0, 1), Point_3(0.5, 0.5, 0.5),
  //                 m);
  // make_tetrahedron(Point_3(1, 0, 0), Point_3(1, 1, 0), Point_3(0, 1, 0), Point_3(0, 0, 0),
  //                 m);
  make_hexahedron(Point_3(1,0,0), Point_3(1,1,0), Point_3(0,1,0), Point_3(0,0,0),
                  Point_3(0,0,1), Point_3(1,0,1), Point_3(1,1,1), Point_3(0,1,1),
                  m);

  timer.start();
  EMesh kernel = PMP::experimental::kernel(m);
  timer.stop();

  std::ofstream("kernel.off") << kernel;
  std::cout << "test_kernel done in " << timer.time() << "\n";
}


void test_kernel(std::string fname)
{
  CGAL::Real_timer timer;
  Mesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }

  timer.start();
  Mesh kernel = PMP::experimental::kernel(m);
  timer.stop();

  std::ofstream("kernel.off") << kernel;
  std::cout << "test_kernel done in " << timer.time() << "\n";
}

void test_kernel_with_rounding(std::string fname)
{
  CGAL::Real_timer timer;
  Mesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }
  round_mesh(m);
  timer.start();
  Mesh kernel = PMP::experimental::kernel(m);
  timer.stop();

  std::ofstream("kernel.off") << kernel;
  std::cout << "test_kernel done in " << timer.time() << "\n";
}

void test_exact_kernel(std::string fname)
{
  CGAL::Real_timer timer;
  EMesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }
  timer.start();
  EMesh kernel = PMP::experimental::kernel(m);
  timer.stop();

  std::ofstream("ekernel.off") << kernel;
  std::cout << "test_exact_kernel done in " << timer.time() << "\n";
}

void test_exact_kernel_with_rounding(std::string fname)
{
  CGAL::Real_timer timer;
  EMesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }
  // round_mesh(m);
  to_integer_mesh(m);
  timer.start();
  EMesh kernel = PMP::experimental::kernel(m, CGAL::parameters::shuffle_planes(true));
  // EMesh kernel = PMP::experimental::kernel(m);
  timer.stop();

  std::ofstream("ekernel.off") << kernel;
  std::cout << "test_exact_kernel done in " << timer.time() << "\n";

  // timer.reset();
  // timer.start();
  // PMP::experimental::kernel(m, CGAL::parameters::shuffle_planes(true));
  // timer.stop();

  // std::cout << "test_exact_kernel with shuffle done in " << timer.time() << "\n";

  // timer.reset();
  // timer.start();
  // PMP::experimental::kernel(m, CGAL::parameters::use_bounding_box_filtering(false));
  // timer.stop();

  // std::cout << "test_exact_kernel without bbox done in " << timer.time() << "\n";

  // timer.reset();
  // timer.start();
  // PMP::experimental::kernel(m, CGAL::parameters::remove_duplicate_planes(true));
  // timer.stop();

  // std::cout << "test_exact_kernel with remove duplicates done in " << timer.time() << "\n";

  // timer.reset();
  // timer.start();
  // PMP::experimental::kernel(m, CGAL::parameters::look_concave_planes_first(true));
  // timer.stop();

  // std::cout << "test_exact_kernel with concave optim done in " << timer.time() << "\n";

  // timer.reset();
  // timer.start();
  // PMP::experimental::kernel(m, CGAL::parameters::look_concave_planes_first(true).remove_duplicate_planes(true));
  // timer.stop();

  // std::cout << "test_exact_kernel with concave optim and plane remove duplicates done in " << timer.time() << "\n";
}

void test_kernel_with_chull(std::string fname)
{
  CGAL::Real_timer timer;
  Mesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }

  timer.start();
  Mesh kernel = PMP::experimental::kernel_using_chull(m);
  timer.stop();

  std::ofstream("kernel_with_chull.off") << kernel;
  std::cout << "test_kernel_with_chull done in " << timer.time() << "\n";
}
void test_kernel_with_chull_and_constructions(std::string fname)
{
  CGAL::Real_timer timer;
  Mesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }

  timer.start();
  Mesh kernel = PMP::experimental::kernel_using_chull_and_constructions(m);
  timer.stop();

  std::ofstream("kernel_with_chull_and_constructions.off") << kernel;
  std::cout << "test_kernel_with_chull_and_constructions done in " << timer.time() << "\n";
}

void test_trettner_kernel(std::string fname)
{
  using int256 = boost::multiprecision::int256_t;

  CGAL::Real_timer timer;
  Mesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }
  to_integer_mesh(m);
  timer.start();
  auto kernel = PMP::experimental::plane_based_kernel(m, CGAL::parameters::geom_traits(CGAL::Homogeneous<int256>()));
  timer.stop();

  std::ofstream("tkernel.off") << kernel;
  std::cout << "test_trettner_kernel done in " << timer.time() << "\n";
}

void test_plane_based_kernel(std::string fname)
{
  CGAL::Real_timer timer;
  EMesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }
  to_integer_mesh(m);
  timer.start();
  auto kernel = PMP::experimental::plane_based_kernel(m);
  timer.stop();

  std::ofstream("tkernel.off") << kernel;
  std::cout << "test plane based kernel with epeck done in " << timer.time() << "\n";
}

int main(int argc, char** argv)
{
  // test_traits();
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  // elementary_test_kernel();
  // test_kernel(filename);
  // test_kernel_with_rounding(filename);
  // test_exact_kernel(filename);
  test_exact_kernel_with_rounding(filename);
  // test_trettner_kernel(filename);
  // test_plane_based_kernel(filename);
  // test_kernel_with_chull(filename);
  // test_kernel_with_chull_and_constructions(filename);
}
