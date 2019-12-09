#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Eigen_linear_algebra_traits.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/optimal_bounding_box.h>
#include <CGAL/Optimal_bounding_box/helper.h>
#include <CGAL/Optimal_bounding_box/nelder_mead_functions.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}

void test_nelder_mead()
{
  typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;
  typedef Linear_algebra_traits::Matrix Matrix;

  Matrix data_points(4,3);
  data_points(0,0) = 0.866802;
  data_points(0,1) = 0.740808;
  data_points(0,2) = 0.895304;
  data_points(1,0) = 0.912651;
  data_points(1,1) = 0.761565;
  data_points(1,2) = 0.160330;
  data_points(2,0) = 0.093661;
  data_points(2,1) = 0.892578;
  data_points(2,2) = 0.737412;
  data_points(3,0) = 0.166461;
  data_points(3,1) = 0.149912;
  data_points(3,2) = 0.364944;

  // one simplex
  std::vector<Matrix> simplex(4);
  Matrix v0(3,3);
  Matrix v1(3,3);
  Matrix v2(3,3);
  Matrix v3(3,3);

  v0(0,0) = -0.2192721;
  v0(0,1) = 0.2792986;
  v0(0,2) = -0.9348326;
  v0(1,0) = -0.7772152;
  v0(1,1) = -0.6292092;
  v0(1,2) = -0.0056861;
  v0(2,0) = -0.5897934;
  v0(2,1) = 0.7253193;
  v0(2,2) = 0.3550431;

  v1(0,0) = -0.588443;
  v1(0,1) = 0.807140;
  v1(0,2) = -0.047542;
  v1(1,0) = -0.786228;
  v1(1,1) = -0.584933;
  v1(1,2) = -0.199246;
  v1(2,0) = -0.188629;
  v1(2,1) = -0.079867;
  v1(2,2) = 0.978795;

  v2(0,0) = -0.277970;
  v2(0,1) = 0.953559;
  v2(0,2) = 0.116010;
  v2(1,0) = -0.567497;
  v2(1,1) = -0.065576;
  v2(1,2) = -0.820760;
  v2(2,0) = -0.775035;
  v2(2,1) = -0.293982;
  v2(2,2) = 0.559370;

  v3(0,0) = -0.32657;
  v3(0,1) = -0.60013;
  v3(0,2) = -0.73020;
  v3(1,0) = -0.20022;
  v3(1,1) = -0.71110;
  v3(1,2) = 0.67398;
  v3(2,0) = -0.92372;
  v3(2,1) = 0.36630;
  v3(2,2) = 0.11207;

  simplex[0] = v0;
  simplex[1] = v1;
  simplex[2] = v2;
  simplex[3] = v3;

  std::size_t nm_iterations = 19;
  CGAL::Optimal_bounding_box::nelder_mead<Linear_algebra_traits>(simplex, data_points, nm_iterations);

  double epsilon = 1e-5;
  const Matrix& v0_new = simplex[0];
  assert(assert_doubles(v0_new(0,0), -0.288975, epsilon));
  assert(assert_doubles(v0_new(0,1), 0.7897657, epsilon));
  assert(assert_doubles(v0_new(0,2), -0.541076, epsilon));
  assert(assert_doubles(v0_new(1,0), -0.9407046, epsilon));
  assert(assert_doubles(v0_new(1,1), -0.3391466, epsilon));
  assert(assert_doubles(v0_new(1,2), 0.0073817, epsilon));
  assert(assert_doubles(v0_new(2,0), -0.1776743, epsilon));
  assert(assert_doubles(v0_new(2,1), 0.5111260, epsilon));
  assert(assert_doubles(v0_new(2,2), 0.84094, epsilon));

  const Matrix& v1_new = simplex[1];
  assert(assert_doubles(v1_new(0,0), -0.458749, epsilon));
  assert(assert_doubles(v1_new(0,1), 0.823283, epsilon));
  assert(assert_doubles(v1_new(0,2), -0.334296, epsilon));
  assert(assert_doubles(v1_new(1,0), -0.885235, epsilon));
  assert(assert_doubles(v1_new(1,1), -0.455997, epsilon));
  assert(assert_doubles(v1_new(1,2), 0.091794, epsilon));
  assert(assert_doubles(v1_new(2,0), -0.076866, epsilon));
  assert(assert_doubles(v1_new(2,1), 0.338040, epsilon));
  assert(assert_doubles(v1_new(2,2), 0.937987, epsilon));

  const Matrix& v2_new = simplex[2];
  assert(assert_doubles(v2_new(0,0), -0.346582, epsilon));
  assert(assert_doubles(v2_new(0,1), 0.878534, epsilon));
  assert(assert_doubles(v2_new(0,2), -0.328724, epsilon));
  assert(assert_doubles(v2_new(1,0), -0.936885, epsilon));
  assert(assert_doubles(v2_new(1,1), -0.341445, epsilon));
  assert(assert_doubles(v2_new(1,2), 0.075251, epsilon));
  assert(assert_doubles(v2_new(2,0), -0.046131, epsilon));
  assert(assert_doubles(v2_new(2,1), 0.334057, epsilon));
  assert(assert_doubles(v2_new(2,2), 0.941423, epsilon));

  const Matrix& v3_new = simplex[3];
  assert(assert_doubles(v3_new(0,0), -0.394713, epsilon));
  assert(assert_doubles(v3_new(0,1), 0.791782, epsilon));
  assert(assert_doubles(v3_new(0,2), -0.466136, epsilon));
  assert(assert_doubles(v3_new(1,0), -0.912112, epsilon));
  assert(assert_doubles(v3_new(1,1), -0.398788, epsilon));
  assert(assert_doubles(v3_new(1,2), 0.094972, epsilon));
  assert(assert_doubles(v3_new(2,0), -0.110692, epsilon));
  assert(assert_doubles(v3_new(2,1), 0.462655, epsilon));
  assert(assert_doubles(v3_new(2,2), 0.879601, epsilon));
}

void test_genetic_algorithm()
{
  CGAL::Eigen_dense_matrix<double, -1, -1> data_points(4, 3); // -1 = dynamic size at run time
  data_points(0,0) = 0.866802;
  data_points(0,1) = 0.740808;
  data_points(0,2) = 0.895304;
  data_points(1,0) = 0.912651;
  data_points(1,1) = 0.761565;
  data_points(1,2) = 0.160330;
  data_points(2,0) = 0.093661;
  data_points(2,1) = 0.892578;
  data_points(2,2) = 0.737412;
  data_points(3,0) = 0.166461;
  data_points(3,1) = 0.149912;
  data_points(3,2) = 0.364944;

  typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;
  CGAL::Optimal_bounding_box::Population<Linear_algebra_traits> pop(5);
  CGAL::Optimal_bounding_box::Evolution<Linear_algebra_traits> evolution(pop, data_points);
  evolution.genetic_algorithm();
  assert(pop.size() == 5);
}

void test_random_unit_tetra()
{
  // this is dynamic at run times
  CGAL::Eigen_dense_matrix<double, -1, -1> data_points(4, 3);

  // points are on their convex hull
  data_points(0,0) = 0.866802;
  data_points(0,1) = 0.740808;
  data_points(0,2) = 0.895304;
  data_points(1,0) = 0.912651;
  data_points(1,1) = 0.761565;
  data_points(1,2) = 0.160330;
  data_points(2,0) = 0.093661;
  data_points(2,1) = 0.892578;
  data_points(2,2) = 0.737412;
  data_points(3,0) = 0.166461;
  data_points(3,1) = 0.149912;
  data_points(3,2) = 0.364944;

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

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  std::ofstream out("data/random_unit_tetra.off");
  out << mesh;
  out.close();
#endif

  typedef CGAL::Eigen_linear_algebra_traits               Linear_algebra_traits;
  typedef Linear_algebra_traits::Matrix                   Matrix;

  std::size_t generations = 10;
  CGAL::Optimal_bounding_box::Population<Linear_algebra_traits> pop(50);
  CGAL::Optimal_bounding_box::Evolution<Linear_algebra_traits> evolution(pop, data_points);
  evolution.evolve(generations);

  Matrix R = evolution.get_best();
  const double epsilon = 1e-3;
  assert(assert_doubles(Linear_algebra_traits::compute_determinant(R), 1, epsilon));
  assert(assert_doubles(R(0,0), -0.25791, epsilon));
  assert(assert_doubles(R(0,1), 0.796512, epsilon));
  assert(assert_doubles(R(0,2), -0.546855, epsilon));
  assert(assert_doubles(R(1,0), -0.947128, epsilon));
  assert(assert_doubles(R(1,1), -0.320242, epsilon));
  assert(assert_doubles(R(1,2), -0.0197553, epsilon));
  assert(assert_doubles(R(2,0), -0.190861, epsilon));
  assert(assert_doubles(R(2,1), 0.512847, epsilon));
  assert(assert_doubles(R(2,2), 0.836992, epsilon));

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  // postprocessing
  CGAL::Eigen_dense_matrix<double> obb(8, 3);
  CGAL::Optimal_bounding_box::post_processing(data_points, R, obb);
#endif
}

void test_reference_tetrahedron(const char* fname)
{
  std::ifstream input(fname);
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  typedef CGAL::Eigen_linear_algebra_traits               Linear_algebra_traits;
  typedef Linear_algebra_traits::Matrix                   Matrix;

  // points in a matrix
  Matrix points;
  CGAL::Optimal_bounding_box::sm_to_matrix(mesh, points);

  std::size_t generations = 10;
  CGAL::Optimal_bounding_box::Population<Linear_algebra_traits> pop(50);
  CGAL::Optimal_bounding_box::Evolution<Linear_algebra_traits> experiment(pop, points);
  experiment.evolve(generations);

  Matrix R = experiment.get_best();
  double epsilon = 1e-5;
  assert(assert_doubles(Linear_algebra_traits::compute_determinant(R), 1, epsilon));

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  // postprocessing
  Matrix obb(8, 3);
  CGAL::Optimal_bounding_box::post_processing(points, R, obb);
#endif
}

void test_long_tetrahedron(const std::string fname)
{
  std::ifstream input(fname);
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  typedef CGAL::Eigen_linear_algebra_traits               Linear_algebra_traits;
  typedef Linear_algebra_traits::Matrix                   Matrix;

  // points in a matrix
  Matrix points;
  CGAL::Optimal_bounding_box::sm_to_matrix(mesh, points);

  std::size_t max_generations = 10;
  CGAL::Optimal_bounding_box::Population<Linear_algebra_traits> pop(50);
  CGAL::Optimal_bounding_box::Evolution<Linear_algebra_traits> experiment(pop, points);
  experiment.evolve(max_generations);

  Matrix R = experiment.get_best();
  double epsilon = 1e-3;
  assert(assert_doubles(Linear_algebra_traits::compute_determinant(R), 1, epsilon));
  assert(assert_doubles(R(0,0), -1, epsilon));
  assert(assert_doubles(R(0,1), 0, epsilon));
  assert(assert_doubles(R(0,2), 0, epsilon));
  assert(assert_doubles(R(1,0), 0, epsilon));
  assert(assert_doubles(R(1,1), -0.707107, epsilon));
  assert(assert_doubles(R(1,2), 0.707106, epsilon) ||
         assert_doubles(R(1,2), -0.707106, epsilon));
  assert(assert_doubles(R(2,0), 0, epsilon));
  assert(assert_doubles(R(2,1), 0.707106, epsilon) ||
         assert_doubles(R(1,2), -0.707106, epsilon));
  assert(assert_doubles(R(2,2), 0.707107, epsilon));

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  // postprocessing
  Matrix obb(8, 3);
  CGAL::Optimal_bounding_box::post_processing(points, R, obb);
#endif
}

void test_compute_obb_evolution(const std::string fname)
{
  std::ifstream input(fname);
  typedef CGAL::Surface_mesh<K::Point_3> SMesh;
  SMesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  // get mesh points
  std::vector<K::Point_3> sm_points;

  typedef typename boost::graph_traits<SMesh>::vertex_descriptor                  vertex_descriptor;
  typedef typename boost::property_map<SMesh, boost::vertex_point_t>::const_type  PointPMap;

  PointPMap pmap = get(boost::vertex_point, mesh);
  for(vertex_descriptor v : vertices(mesh))
    sm_points.push_back(get(pmap, v));


  CGAL::Eigen_linear_algebra_traits la_traits;
  std::vector<K::Point_3> obb_points;
  CGAL::Optimal_bounding_box::compute_optimal_bounding_box(sm_points, obb_points, la_traits, true);

  double epsilon = 1e-3;
  double vol = CGAL::Optimal_bounding_box::calculate_volume(obb_points);
  assert(assert_doubles(vol, 0.883371, epsilon));

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  /*
  for(int i = 0; i < 8; ++i)
    std::cout << obb_points[i].x() << " " << obb_points[i].y() << " " << obb_points[i].z() << "\n" ;
  */
  CGAL::Surface_mesh<K::Point_3> result_mesh;
  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
      obb_points[4], obb_points[5], obb_points[6], obb_points[7], result_mesh);

  std::ofstream out("data/obb_result.off");
  out << result_mesh;
  out.close();
#endif
}

void test_compute_obb_mesh(const std::string fname)
{
  std::ifstream input(fname);
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  CGAL::Eigen_linear_algebra_traits la_traits;
  CGAL::Surface_mesh< K::Point_3> obbmesh;
  CGAL::Optimal_bounding_box::compute_optimal_bounding_box(mesh, obbmesh, la_traits, true);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  std::ofstream out("/tmp/result_elephant.off");
  out << obbmesh;
  out.close();
#endif
}

void test_function_defaults_traits(const std::string fname1,
                                   const std::string fname2)
{
  std::ifstream input1(fname1);
  CGAL::Surface_mesh<K::Point_3> mesh1;
  if (!input1 || !(input1 >> mesh1) || mesh1.is_empty())
  {
    std::cerr << fname1 << " is not a valid off file.\n";
    std::exit(1);
  }

  std::ifstream input2(fname2);
  CGAL::Surface_mesh<K::Point_3> mesh2;
  if (!input2 || !(input2 >> mesh2) || mesh2.is_empty())
  {
    std::cerr << fname2 << " is not a valid off file.\n";
    std::exit(1);
  }

  // test one
  std::vector<K::Point_3> sm_points;
  typedef CGAL::Surface_mesh<K::Point_3> SMesh;

  typedef typename boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_map<SMesh, boost::vertex_point_t>::const_type PointPMap;
  PointPMap pmap = get(boost::vertex_point, mesh1);

  for(vertex_descriptor v : vertices(mesh1))
    sm_points.push_back(get(pmap, v));

  std::vector<K::Point_3> obb_points;
  CGAL::Optimal_bounding_box::compute_optimal_bounding_box(sm_points, obb_points, true);

  const double epsilon = 1e-3;
  const double vol = CGAL::Optimal_bounding_box::calculate_volume(obb_points);
  assert(assert_doubles(vol, 0.883371, epsilon));

  // test two
  CGAL::Surface_mesh<K::Point_3> obbmesh;
  CGAL::Optimal_bounding_box::compute_optimal_bounding_box(mesh2, obbmesh, true);
}

int main()
{
  test_nelder_mead();
  test_genetic_algorithm();
  test_random_unit_tetra();
  test_reference_tetrahedron("data/reference_tetrahedron.off");
  test_long_tetrahedron("data/long_tetrahedron.off");
  test_compute_obb_evolution("data/random_unit_tetra.off");
  test_compute_obb_mesh("data/elephant.off");
  test_function_defaults_traits("data/random_unit_tetra.off", "data/elephant.off");

  return 0;
}
