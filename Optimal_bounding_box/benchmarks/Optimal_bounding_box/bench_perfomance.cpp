#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Eigen_linear_algebra_traits.h>
#include <CGAL/Optimal_bounding_box/optimal_bounding_box.h>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Timer.h>

#include <fstream>
#include <iostream>

//#define CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

void bench_finding_obb(std::string fname)
{
  std::ifstream input(fname);

  // import a mesh
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  // export some times
  std::ofstream outt("data/times.txt");
  outt << "nb_vertices "<<  "with_ch " << "without_ch" << std::endl;

  CGAL::Timer timer;
  std::size_t measurements = 4;
  CGAL::Eigen_linear_algebra_traits la_traits;

  for (std::size_t t = 0; t < measurements; ++t)
  {
    std::cout << "#vertices= " << vertices(mesh).size() << " |";

    // 1) using convex hull
    timer.start();
    CGAL::Surface_mesh<K::Point_3> obb_mesh1;
    CGAL::Optimal_bounding_box::compute_optimal_bounding_box(mesh, obb_mesh1, la_traits, true);
    timer.stop();
    double t_ch = timer.time();
    std::cout << " with ch: " << timer.time() << " s |";

    // 2) without convex hull
    timer.reset();
    timer.start();
    CGAL::Surface_mesh<K::Point_3> obb_mesh2;
    CGAL::Optimal_bounding_box::compute_optimal_bounding_box(mesh, obb_mesh2, la_traits, false);
    timer.stop();
    double t_no_ch = timer.time();
    std::cout << " without ch: " <<  timer.time() << " s\n";
    timer.reset();

    outt << vertices(mesh).size() << " " << t_ch << " " << t_no_ch << std::endl;

    // 3) subdivision
    CGAL::Subdivision_method_3::CatmullClark_subdivision(mesh,
                                      CGAL::parameters::number_of_iterations(1));
  }

  outt.close();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK
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
