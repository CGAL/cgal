#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Optimal_bounding_box/optimization_algorithms.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/obb.h>
#include <iostream>
#include <fstream>

#include <CGAL/Timer.h>

//#define OBB_DEBUG_BENCHMARK


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}

/*
template <typename SurfaceMesh, typename Matrix>
void sm_to_matrix(SurfaceMesh& sm, Matrix& mat)
{
  typedef typename boost::property_map<SurfaceMesh, boost::vertex_point_t>::const_type Vpm;
  typedef typename boost::property_traits<Vpm>::reference Point_ref;
  typedef typename boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;
  Vpm vpm = get(boost::vertex_point, sm);

  mat.resize(vertices(sm).size(), 3);
  std::size_t i = 0;
  for(vertex_descriptor v : vertices(sm))
  {
    Point_ref p = get(vpm, v);
    mat(i, 0) = p.x();
    mat(i, 1) = p.y();
    mat(i, 2) = p.z();
    ++i;
  }
}
*/

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

  std::vector<K::Point_3> sm_points;

  // import a mesh
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }
  // get mesh points
  gather_mesh_points(mesh, sm_points);

  CGAL::Timer timer;
  timer.start();

  // 1) using convex hull
  std::vector<K::Point_3> obb_points1;
  CGAL::Optimal_bounding_box::find_obb(sm_points, obb_points1, true);

  timer.stop();
  std::cout << "found obb with convex hull: " << timer.time() << " seconds\n";

  timer.reset();
  timer.start();

  // 2) without convex hull
  std::vector<K::Point_3> obb_points2;
  CGAL::Optimal_bounding_box::find_obb(sm_points, obb_points2, false);

  timer.stop();

  std::cout << "found obb without convex hull: " <<  timer.time() << " seconds\n";
  timer.reset();


  double epsilon = 1e-3;
  double vol1 = calculate_volume(obb_points1);
  double vol2 = calculate_volume(obb_points2);
  //std::cout << "vol1= " << vol1 << " -- " << "vol2= " << vol2 << std::endl;


#ifdef OBB_DEBUG_BENCHMARK
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

int main(int argc, char* argv[])
{
  bench_finding_obb("data/elephant.off");

  return 0;
}
