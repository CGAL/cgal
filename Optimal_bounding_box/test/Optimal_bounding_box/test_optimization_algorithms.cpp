#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Optimal_bounding_box/optimization_algorithms.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/obb.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::Matrix3d Matrix3d;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}

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
    mat.coeffRef(i, 0) = p.x();
    mat.coeffRef(i, 1) = p.y();
    mat.coeffRef(i, 2) = p.z();
    ++i;
  }
}

// it is called after post processing
template <typename Matrix>
void matrix_to_mesh_and_draw(Matrix& data_points, std::string filename)
{
  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> Mesh;

  // Simplex -> std::vector
  std::vector<Point> points;

  for(int i = 0; i < data_points.rows(); ++i)
  {
    Point p(data_points(i, 0), data_points(i, 1), data_points(i, 2));
    points.push_back(p);
  }

  Mesh mesh;
  CGAL::make_hexahedron(points[0], points[1], points[2], points[3], points[4], points[5],
      points[6], points[7], mesh);

  std::ofstream out(filename);
  out << mesh;
  out.close();
}

template <typename Point>
double calculate_volume(std::vector<Point> points)
{
  CGAL::Bbox_3 bbox;
  bbox = bbox_3(points.begin(), points.end());
  K::Iso_cuboid_3 ic(bbox);
  return ic.volume();
}

void test_population()
{
  CGAL::Optimal_bounding_box::Population<Matrix3d> pop(5);
  //pop.show_population();
  CGAL_assertion(pop.size() == 5);
}

void test_nelder_mead()
{
  Eigen::Matrix<double, 4, 3> data_points(4, 3);
  data_points << 0.866802, 0.740808, 0.895304,
                 0.912651, 0.761565, 0.160330,
                 0.093661, 0.892578, 0.737412,
                 0.166461, 0.149912, 0.364944;

  // one simplex
  std::vector<Matrix3d> simplex(4);

  Matrix3d v0;
  v0 <<  -0.2192721,   0.2792986,  -0.9348326,
          -0.7772152,  -0.6292092,  -0.0056861,
          -0.5897934,   0.7253193,   0.3550431;

  Matrix3d v1;
  v1 <<  -0.588443,   0.807140,  -0.047542,
          -0.786228,  -0.584933,  -0.199246,
          -0.188629,  -0.079867,   0.978795;

  Matrix3d v2;
  v2 << -0.277970,  0.953559,  0.116010,
         -0.567497,  -0.065576,   -0.820760,
         -0.775035,   -0.293982,   0.559370;

  Matrix3d v3;
  v3 <<   -0.32657,  -0.60013,  -0.73020,
           -0.20022,  -0.71110,   0.67398,
           -0.92372,   0.36630,   0.11207;

  simplex[0] = v0;
  simplex[1] = v1;
  simplex[2] = v2;
  simplex[3] = v3;

  std::size_t iterations = 19;
  CGAL::Optimal_bounding_box::nelder_mead(simplex, data_points, iterations);

  double epsilon = 1e-5;

  Matrix3d v0_new = simplex[0];
  CGAL_assertion(assert_doubles(v0_new(0,0), -0.288975, epsilon));
  CGAL_assertion(assert_doubles(v0_new(0,1), 0.7897657, epsilon));
  CGAL_assertion(assert_doubles(v0_new(0,2), -0.541076, epsilon));
  CGAL_assertion(assert_doubles(v0_new(1,0), -0.9407046, epsilon));
  CGAL_assertion(assert_doubles(v0_new(1,1), -0.3391466, epsilon));
  CGAL_assertion(assert_doubles(v0_new(1,2), 0.0073817, epsilon));
  CGAL_assertion(assert_doubles(v0_new(2,0), -0.1776743, epsilon));
  CGAL_assertion(assert_doubles(v0_new(2,1), 0.5111260, epsilon));
  CGAL_assertion(assert_doubles(v0_new(2,2), 0.84094, epsilon));

  Matrix3d v1_new = simplex[1];
  CGAL_assertion(assert_doubles(v1_new(0,0), -0.458749, epsilon));
  CGAL_assertion(assert_doubles(v1_new(0,1), 0.823283, epsilon));
  CGAL_assertion(assert_doubles(v1_new(0,2), -0.334296, epsilon));
  CGAL_assertion(assert_doubles(v1_new(1,0), -0.885235, epsilon));
  CGAL_assertion(assert_doubles(v1_new(1,1), -0.455997, epsilon));
  CGAL_assertion(assert_doubles(v1_new(1,2), 0.091794, epsilon));
  CGAL_assertion(assert_doubles(v1_new(2,0), -0.076866, epsilon));
  CGAL_assertion(assert_doubles(v1_new(2,1), 0.338040, epsilon));
  CGAL_assertion(assert_doubles(v1_new(2,2), 0.937987, epsilon));

  Matrix3d v2_new = simplex[2];
  CGAL_assertion(assert_doubles(v2_new(0,0), -0.346582, epsilon));
  CGAL_assertion(assert_doubles(v2_new(0,1), 0.878534, epsilon));
  CGAL_assertion(assert_doubles(v2_new(0,2), -0.328724, epsilon));
  CGAL_assertion(assert_doubles(v2_new(1,0), -0.936885, epsilon));
  CGAL_assertion(assert_doubles(v2_new(1,1), -0.341445, epsilon));
  CGAL_assertion(assert_doubles(v2_new(1,2), 0.075251, epsilon));
  CGAL_assertion(assert_doubles(v2_new(2,0), -0.046131, epsilon));
  CGAL_assertion(assert_doubles(v2_new(2,1), 0.334057, epsilon));
  CGAL_assertion(assert_doubles(v2_new(2,2), 0.941423, epsilon));

  Matrix3d v3_new = simplex[3];
  CGAL_assertion(assert_doubles(v3_new(0,0), -0.394713, epsilon));
  CGAL_assertion(assert_doubles(v3_new(0,1), 0.791782, epsilon));
  CGAL_assertion(assert_doubles(v3_new(0,2), -0.466136, epsilon));
  CGAL_assertion(assert_doubles(v3_new(1,0), -0.912112, epsilon));
  CGAL_assertion(assert_doubles(v3_new(1,1), -0.398788, epsilon));
  CGAL_assertion(assert_doubles(v3_new(1,2), 0.094972, epsilon));
  CGAL_assertion(assert_doubles(v3_new(2,0), -0.110692, epsilon));
  CGAL_assertion(assert_doubles(v3_new(2,1), 0.462655, epsilon));
  CGAL_assertion(assert_doubles(v3_new(2,2), 0.879601, epsilon));
}

void test_genetic_algorithm()
{
  Eigen::Matrix<double, 4, 3> data_points(4, 3);
  data_points << 0.866802, 0.740808, 0.895304,
                 0.912651, 0.761565, 0.160330,
                 0.093661, 0.892578, 0.737412,
                 0.166461, 0.149912, 0.364944;

  CGAL::Optimal_bounding_box::Population<Matrix3d> pop(5);
  CGAL::Optimal_bounding_box::genetic_algorithm(pop, data_points);
  CGAL_assertion(pop.size() == 5);
}

void test_random_unit_tetra()
{
  Eigen::Matrix<double, 4, 3> data_points(4, 3); // points on their convex hull
  data_points << 0.866802, 0.740808, 0.895304,
                 0.912651, 0.761565, 0.160330,
                 0.093661, 0.892578, 0.737412,
                 0.166461, 0.149912, 0.364944;

  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> Mesh;

  // make a mesh and export it
  Mesh mesh;
  Point p1(0.866802, 0.740808, 0.895304);
  Point p2(0.912651, 0.761565, 0.160330);
  Point p3(0.093661, 0.892578, 0.737412);
  Point p4(0.166461, 0.149912, 0.364944);
  CGAL::make_tetrahedron(p1, p2, p3, p4, mesh);

#ifdef OBB_DEBUG_TEST
  std::ofstream out("data/random_unit_tetra.off");
  out << mesh;
  out.close();
#endif

  Matrix3d R;
  std::size_t generations = 10;
  CGAL::Optimal_bounding_box::evolution(R, data_points, generations);

  double epsilon = 1e-3;
  CGAL_assertion(assert_doubles(R.determinant(), 1, epsilon));
  CGAL_assertion(assert_doubles(R(0,0), -0.25791, epsilon));
  CGAL_assertion(assert_doubles(R(0,1), 0.796512, epsilon));
  CGAL_assertion(assert_doubles(R(0,2), -0.546855, epsilon));
  CGAL_assertion(assert_doubles(R(1,0), -0.947128, epsilon));
  CGAL_assertion(assert_doubles(R(1,1), -0.320242, epsilon));
  CGAL_assertion(assert_doubles(R(1,2), -0.0197553, epsilon));
  CGAL_assertion(assert_doubles(R(2,0), -0.190861, epsilon));
  CGAL_assertion(assert_doubles(R(2,1), 0.512847, epsilon));
  CGAL_assertion(assert_doubles(R(2,2), 0.836992, epsilon));

#ifdef OBB_DEBUG_TEST
  // postprocessing
  MatrixXd obb(8, 3);
  CGAL::Optimal_bounding_box::post_processing(data_points, R, obb);
  matrix_to_mesh_and_draw(obb, "data/random_unit_tetra_result.off");
#endif
}

void test_reference_tetrahedron(const char* fname)
{
  std::ifstream input(fname);
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  // points in a matrix
  MatrixXd points;
  sm_to_matrix(mesh, points);

  Matrix3d R(3, 3);
  std::size_t generations = 10;
  CGAL::Optimal_bounding_box::evolution(R, points, generations);
  double epsilon = 1e-5;
  CGAL_assertion(assert_doubles(R.determinant(), 1, epsilon));

  #ifdef OBB_DEBUG_TEST
  // postprocessing
  MatrixXd obb(8, 3);
  CGAL::Optimal_bounding_box::post_processing(points, R, obb);
  matrix_to_mesh_and_draw(obb, "data/OBB.off");
  #endif
}

void test_long_tetrahedron(std::string fname)
{
  std::ifstream input(fname);
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  // points in a matrix
  MatrixXd points;
  sm_to_matrix(mesh, points);

  MatrixXd R(3, 3); // test with dynamic size
  std::size_t generations = 10;
  CGAL::Optimal_bounding_box::evolution(R, points, generations);
  double epsilon = 1e-3;
  CGAL_assertion(assert_doubles(R.determinant(), 1, epsilon));
  CGAL_assertion(assert_doubles(R(0,0), -1, epsilon));
  CGAL_assertion(assert_doubles(R(0,1), 0, epsilon));
  CGAL_assertion(assert_doubles(R(0,2), 0, epsilon));
  CGAL_assertion(assert_doubles(R(1,0), 0, epsilon));
  CGAL_assertion(assert_doubles(R(1,1), -0.707107, epsilon));
  CGAL_assertion(assert_doubles(R(1,2), 0.707106, epsilon) ||
                 assert_doubles(R(1,2), -0.707106, epsilon));
  CGAL_assertion(assert_doubles(R(2,0), 0, epsilon));
  CGAL_assertion(assert_doubles(R(2,1), 0.707106, epsilon) ||
                 assert_doubles(R(1,2), -0.707106, epsilon));
  CGAL_assertion(assert_doubles(R(2,2), 0.707107, epsilon));

  #ifdef OBB_DEBUG_TEST
  // postprocessing
  MatrixXd obb(8, 3);
  CGAL::Optimal_bounding_box::post_processing(points, R, obb);
  matrix_to_mesh_and_draw(obb, fname + "result.off");
  #endif
}

void test_find_obb(std::string fname)
{
  std::ifstream input(fname);
  typedef CGAL::Surface_mesh<K::Point_3> SMesh;
  SMesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  // get mesh points
  std::vector<K::Point_3> sm_points;
  typedef typename boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_map<SMesh, boost::vertex_point_t>::const_type PointPMap;
  PointPMap pmap = get(boost::vertex_point, mesh);
  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
    sm_points.push_back(get(pmap, v));

  std::vector<K::Point_3> obb_points;
  CGAL::Optimal_bounding_box::find_obb(sm_points, obb_points, true);

  double epsilon = 1e-3;
  double vol = calculate_volume(obb_points);
  CGAL_assertion(assert_doubles(vol, 0.883371, epsilon));

  #ifdef OBB_DEBUG_TEST
  for(int i = 0; i < 8; ++i)
  {
    std::cout << obb_points[i].x() << " " << obb_points[i].y() << " " << obb_points[i].z() << "\n" ;
  }
  CGAL::Surface_mesh<K::Point_3> result_mesh;
  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                        obb_points[4], obb_points[5], obb_points[6], obb_points[7], result_mesh);

  std::ofstream out("data/obb_result.off");
  out << result_mesh;
  out.close();
  #endif
}

void test_find_obb_mesh(std::string fname)
{
  std::ifstream input(fname);
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  CGAL::Surface_mesh<K::Point_3> obbmesh;
  CGAL::Optimal_bounding_box::find_obb(mesh, obbmesh, true);

  #ifdef OBB_DEBUG_TEST
  std::ofstream out("/tmp/result_elephant.off");
  out << obbmesh;
  out.close();
  #endif
}

int main(int argc, char* argv[])
{
  test_population();
  test_nelder_mead();
  test_genetic_algorithm();
  test_random_unit_tetra();
  test_reference_tetrahedron("data/reference_tetrahedron.off");
  test_long_tetrahedron("data/long_tetrahedron.off");
  test_find_obb("data/random_unit_tetra.off");
  test_find_obb_mesh("data/elephant.off");

  return 0;
}
