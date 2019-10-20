#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Eigen_linear_algebra_traits.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/optimal_bounding_box.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

//#define CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}

template <typename SurfaceMesh, typename Point>
void gather_mesh_points(SurfaceMesh& mesh, std::vector<Point>& points)
{
  typedef typename boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_map<SurfaceMesh, CGAL::vertex_point_t>::type PointPMap;
  PointPMap pmap = get(boost::vertex_point, mesh);
  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
    points.push_back(get(pmap, v));
}

template <typename Point>
double calculate_volume(std::vector<Point> points)
{
  CGAL::Bbox_3 bbox;
  bbox = bbox_3(points.begin(), points.end());
  K::Iso_cuboid_3 ic(bbox);
  return ic.volume();
}

void bench_finding_obb(std::string fname)
{
  std::ifstream input(fname);

  CGAL::Eigen_linear_algebra_traits la_traits;
  std::vector<K::Point_3> sm_points;

  // import a mesh
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  // get mesh points
  gather_mesh_points(mesh, sm_points);

  CGAL::Timer timer;

  // 1) measure convex hull calculation
  timer.start();
  CGAL::Polyhedron_3<K> poly;
  convex_hull_3(sm_points.begin(), sm_points.end(), poly);
  std::vector<K::Point_3> ch_points(poly.points_begin(), poly.points_end());
  timer.stop();
  std::cout << "takes : " << timer.time() << " seconds to find the convex hull\n";

  // 2) using convex hull
  timer.reset();
  timer.start();
  std::vector<K::Point_3> obb_points1;
  CGAL::Optimal_bounding_box::compute_optimal_bounding_box(sm_points, obb_points1, la_traits, true);
  timer.stop();
  std::cout << "found obb using convex hull: " << timer.time() << " seconds\n";

  // 3) without convex hull
  timer.reset();
  timer.start();
  std::vector<K::Point_3> obb_points2;
  CGAL::Optimal_bounding_box::compute_optimal_bounding_box(sm_points, obb_points2, la_traits, false);
  timer.stop();
  std::cout << "found obb without convex hull: " <<  timer.time() << " seconds\n";
  timer.reset();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK
  CGAL::Surface_mesh<K::Point_3> result_mesh1;
  CGAL::make_hexahedron(obb_points1[0], obb_points1[1], obb_points1[2], obb_points1[3],
                        obb_points1[4], obb_points1[5], obb_points1[6], obb_points1[7],
                        result_mesh1);

  CGAL::Surface_mesh<K::Point_3> result_mesh2;
  CGAL::make_hexahedron(obb_points2[0], obb_points2[1], obb_points2[2], obb_points2[3],
                        obb_points2[4], obb_points2[5], obb_points2[6], obb_points2[7],
                        result_mesh2);

  std::ofstream out1("data/obb_result1.off");
  out1 << result_mesh1;
  out1.close();

  std::ofstream out2("data/obb_result2.off");
  out2 << result_mesh2;
  out2.close();
#endif
}

int main()
{
  bench_finding_obb("data/elephant.off");

  return 0;
}
