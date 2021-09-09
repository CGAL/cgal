#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Optimal_bounding_box/oriented_bounding_box.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK

typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
typedef K::Point_3                                                     Point_3;

typedef CGAL::Surface_mesh<Point_3>                                    Surface_mesh;
typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor  vertex_descriptor;

void bench_finding_obb(const std::string filename,
                       const int iter)
{
  CGAL::Timer timer;

  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > unused_faces;

  CGAL::IO::read_polygon_soup(filename, points, unused_faces);

  std::vector<Point_3> ch_points;
  std::array<Point_3, 8> obb_points1;
  std::array<Point_3, 8> obb_points2;

  double ch_time = 0.;
  double obb_time = 0.;
  double total_time_ch = 0.;
  double total_time_no_ch = 0.;

  for(int i=0; i<iter; ++i)
  {
//    std::cout << "Iter #" << i << std::endl;
//    std::cout << points.size() << std::endl;

    // 1) measure convex hull calculation
    timer.reset();
    ch_points.clear();

    timer.start();
    extreme_points_3(points, std::back_inserter(ch_points));
    timer.stop();
    ch_time += timer.time();

    // 2) using convex hull
    timer.reset();

    timer.start();
    CGAL::oriented_bounding_box(ch_points, obb_points1, CGAL::parameters::use_convex_hull(false));
    timer.stop();
    obb_time += timer.time();

    // 2bis) using convex hull in one go
    timer.reset();

    timer.start();
    CGAL::oriented_bounding_box(points, obb_points1, CGAL::parameters::use_convex_hull(true));
    timer.stop();
    total_time_ch += timer.time();

    // 3) without convex hull
    timer.reset();

    timer.start();
    CGAL::oriented_bounding_box(points, obb_points2, CGAL::parameters::use_convex_hull(false));
    timer.stop();
    total_time_no_ch += timer.time();

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
    oss1 << "data/obb_result1_iter_" << i << std::ends;
    std::ofstream out1(oss1.str().c_str());
    out1 << result_mesh1;

    std::stringstream oss2;
    oss2 << "data/obb_result2_iter_" << i << std::ends;
    std::ofstream out2(oss2.str().c_str());
    out2 << result_mesh2;
#endif
  }

#if 1 // only outputs the core stuff
  std::cout << points.size() << " "
            << ch_points.size() << " "
            << ch_time / iter << " "
            << obb_time / iter << " "
            << total_time_ch / iter << " "
            << total_time_no_ch / iter << std::endl;
#else
  std::cout << "Average of " << iter << " iterations" << std::endl;
  std::cout << points.size() << " vertices in the mesh" << std::endl;
  std::cout << ch_points.size() << " vertices on the convex hull" << std::endl;
  std::cout << ch_time / iter << " seconds to find the convex hull" << std::endl;
  std::cout << obb_time / iter << " seconds to find the best rotation" << std::endl;
  std::cout << total_time_ch / iter << " seconds to compute and find. "
            << "Should be about equal to: " << (ch_time + obb_time) / iter << std::endl;
  std::cout << total_time_no_ch / iter << " seconds to find the best rotation (NO CH)" << std::endl;
#endif
}

int main(int argc, char** argv)
{
  const int iter = (argc > 2) ? std::atoi(argv[2]) : 5;
  bench_finding_obb((argc > 1) ? argv[1] : "data/elephant.off", iter);

  return EXIT_SUCCESS;
}
