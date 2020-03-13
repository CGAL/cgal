#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Optimal_bounding_box/oriented_bounding_box.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/subdivision_method_3.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK

typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
typedef K::Point_3                                                     Point_3;

typedef CGAL::Surface_mesh<Point_3>                                    Surface_mesh;
typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor  vertex_descriptor;

void bench_finding_obb(const std::string fname)
{
  std::ifstream input(fname);

  // import a mesh
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  CGAL::Timer timer;

  std::size_t measurements = 4;
  for(std::size_t t=0; t<measurements; ++t)
  {
    std::cout << "Iteration: " << t << std::endl;
    std::cout << num_vertices(mesh) << " nv and " << num_faces(mesh) << " nf" << std::endl;

    // 1) measure convex hull calculation
    timer.start();
    CGAL::Polyhedron_3<K> poly;
    convex_hull_3(mesh, poly);
    timer.stop();
    std::cout << "takes : " << timer.time() << " seconds to find the convex hull\n";
    std::cout << num_vertices(poly) << " vertices on the convex hull" << std::endl;

    // 2) using convex hull
    timer.reset();
    timer.start();
    std::array<Point_3, 8> obb_points1;
    CGAL::oriented_bounding_box(mesh, obb_points1, CGAL::parameters::use_convex_hull(true));
    timer.stop();
    std::cout << "found obb using convex hull: " << timer.time() << " seconds\n";

    // 3) without convex hull
    timer.reset();
    timer.start();
    std::array<Point_3, 8> obb_points2;
    CGAL::oriented_bounding_box(mesh, obb_points2, CGAL::parameters::use_convex_hull(false));
    timer.stop();
    std::cout << "found obb without convex hull: " <<  timer.time() << " seconds\n";
    timer.reset();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK
    Surface_mesh result_mesh1;
    CGAL::make_hexahedron(obb_points1[0], obb_points1[1], obb_points1[2], obb_points1[3],
                          obb_points1[4], obb_points1[5], obb_points1[6], obb_points1[7],
                          result_mesh1);

    Surface_mesh result_mesh2;
    CGAL::make_hexahedron(obb_points2[0], obb_points2[1], obb_points2[2], obb_points2[3],
                          obb_points2[4], obb_points2[5], obb_points2[6], obb_points2[7],
                          result_mesh2);

    std::stringstream oss1;
    oss1 << "data/obb_result1_iter_" << t << std::ends;
    std::ofstream out1(oss1.str().c_str());
    out1 << result_mesh1;

    std::stringstream oss2;
    oss2 << "data/obb_result2_iter_" << t << std::ends;
    std::ofstream out2(oss2.str().c_str());
    out2 << result_mesh2;

    CGAL::Subdivision_method_3::CatmullClark_subdivision(mesh, CGAL::parameters::number_of_iterations(1));
#endif
  }
}

int main(int argc, char** argv)
{
  bench_finding_obb((argc > 1) ? argv[1] : "data/elephant.off");

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
