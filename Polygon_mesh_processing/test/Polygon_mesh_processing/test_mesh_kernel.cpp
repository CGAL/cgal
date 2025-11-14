#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/internal/kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Real_timer.h>

using K=CGAL::Exact_predicates_inexact_constructions_kernel;
using EK=CGAL::Exact_predicates_exact_constructions_kernel;
using Mesh=CGAL::Surface_mesh<K::Point_3>;
using EMesh=CGAL::Surface_mesh<EK::Point_3>;

namespace PMP=CGAL::Polygon_mesh_processing;

void test_traits()
{
  Mesh m;
  std::ifstream(CGAL::data_file_path("meshes/elephant.off")) >> m;

  std::pair p(0,0.05);
  using Traits = PMP::Orthogonal_cut_plane_traits<K>;

  PMP::clip(m, p, CGAL::parameters::geom_traits(Traits()));

  std::ofstream("clipped.off") << m;
}

void test_kernel(std::string fname)
{
  CGAL::Real_timer timer;
  Mesh m;
  std::ifstream(fname) >> m;

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
  std::ifstream(fname) >> m;
  timer.start();
  EMesh kernel = PMP::experimental::kernel(m);
  timer.stop();

  std::ofstream("ekernel.off") << kernel;
  std::cout << "test_exact_kernel done in " << timer.time() << "\n";
}

void test_kernel_with_chull(std::string fname)
{
  CGAL::Real_timer timer;
  Mesh m;
  std::ifstream(fname) >> m;

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
  std::ifstream(fname) >> m;

  timer.start();
  Mesh kernel = PMP::experimental::kernel_using_chull_and_constructions(m);
  timer.stop();

  std::ofstream("kernel_with_chull_and_constructions.off") << kernel;
  std::cout << "test_kernel_with_chull_and_constructions done in " << timer.time() << "\n";
}

int main(int argc, char** argv)
{
  test_traits();
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  test_kernel(filename);
  test_exact_kernel(filename);
  test_kernel_with_chull(filename);
  test_kernel_with_chull_and_constructions(filename);
}
