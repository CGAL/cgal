#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/obb.h>
#include <CGAL/Eigen_linear_algebra_traits.h>
#include <iostream>
#include <fstream>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Timer.h>

//#define OBB_DEBUG_BENCHMARK

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

void bench_finding_obb(std::string fname)
{
  std::ifstream input(fname);

  // import a mesh
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  CGAL::Timer timer;
  std::size_t measurements = 4;
  CGAL::Eigen_linear_algebra_traits la_traits;

  for (std::size_t t = 0; t < measurements; ++t)
  {
    std::cout << "#vertices= " << vertices(mesh).size() << " |";

    timer.start();

    // 1) using convex hull
    CGAL::Surface_mesh<K::Point_3> obb_mesh1;
    CGAL::Optimal_bounding_box::find_obb(mesh, obb_mesh1, la_traits, true);

    timer.stop();
    std::cout << " with ch: " << timer.time() << " s |";

    timer.reset();
    timer.start();

    // 2) without convex hull
    CGAL::Surface_mesh<K::Point_3> obb_mesh2;
    CGAL::Optimal_bounding_box::find_obb(mesh, obb_mesh2, la_traits, false);

    timer.stop();

    std::cout << " without ch: " <<  timer.time() << " s\n";
    timer.reset();

    // 3) subdivision
    CGAL::Subdivision_method_3::CatmullClark_subdivision(mesh,
                                      CGAL::parameters::number_of_iterations(1));
  }


#ifdef OBB_DEBUG_BENCHMARK
  std::ofstream out1("data/obb_result1.off");
  out1 << obb_points1;
  out1.close();

  std::ofstream out2("data/obb_result2.off");
  out2 << obb_points2;
  out2.close();
#endif
}

int main()
{
  bench_finding_obb("data/elephant.off");

  return 0;
}
