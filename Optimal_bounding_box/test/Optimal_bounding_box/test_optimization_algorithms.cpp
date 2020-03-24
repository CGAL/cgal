#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/optimal_bounding_box.h>

#include <array>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point_3;

typedef CGAL::Surface_mesh<Point>                                   Mesh;

typedef CGAL::Oriented_bounding_box_traits_3<K>                     Traits;
typedef Traits::Matrix                                              Matrix;

void check_equality(const FT d1, const FT d2)
{
  const FT epsilon = 1e-3;

  bool ok;
  if(std::is_floating_point<FT>::value)
    ok = CGAL::abs(d1 - d2) < epsilon * CGAL::abs(d2);
  else
    ok = (d1 == d2);

  if(!ok)
  {
    std::cout << "Error: got " << d1 << " but expected: " << d2 << std::endl;
    assert(false);
  }
}

void test_genetic_algorithm()
{
  std::array<Point_3, 4> points;
  points[0] = Point_3(0.866802, 0.740808, 0.895304);
  points[1] = Point_3(0.912651, 0.761565, 0.160330);
  points[2] = Point_3(0.093661, 0.892578, 0.737412);
  points[3] = Point_3(0.166461, 0.149912, 0.364944);

  CGAL::Optimal_bounding_box::internal::Population<Traits> pop(5);
  CGAL::Optimal_bounding_box::internal::Evolution<Traits> evolution(pop, data_points);
  evolution.genetic_algorithm();
  assert(pop.size() == 5);
}

void test_random_unit_tetra()
{
  std::array<Point_3, 4> points;
  points[0] = Point_3(0.866802, 0.740808, 0.895304);
  points[1] = Point_3(0.912651, 0.761565, 0.160330);
  points[2] = Point_3(0.093661, 0.892578, 0.737412);
  points[3] = Point_3(0.166461, 0.149912, 0.364944);

  Mesh mesh;
  CGAL::make_tetrahedron(points[0], points[1], points[2], points[3], mesh);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  std::ofstream out("data/random_unit_tetra.off");
  out << mesh;
  out.close();
#endif

  std::size_t generations = 10;
  CGAL::Optimal_bounding_box::internal::Population<Traits> pop(50);
  CGAL::Optimal_bounding_box::internal::Evolution<Traits> evolution(pop, data_points);
  evolution.evolve(generations);

  Matrix R = evolution.get_best();
  check_equality(Traits::compute_determinant(R), 1);
  check_equality(R(0,0), -0.25791);
  check_equality(R(0,1), 0.796512);
  check_equality(R(0,2), -0.546855);
  check_equality(R(1,0), -0.947128);
  check_equality(R(1,1), -0.320242);
  check_equality(R(1,2), -0.0197553);
  check_equality(R(2,0), -0.190861);
  check_equality(R(2,1), 0.512847);
  check_equality(R(2,2), 0.836992);
}

void test_reference_tetrahedron(const char* fname)
{
  std::ifstream input(fname);
  CGAL::Surface_mesh<Point_3> mesh;
  if(!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file." << std::endl;
    std::exit(1);
  }

  // points in a matrix
  Matrix points;
  CGAL::Optimal_bounding_box::sm_to_matrix(mesh, points);

  std::size_t generations = 10;
  CGAL::Optimal_bounding_box::internal::Population<Traits> pop(50);
  CGAL::Optimal_bounding_box::internal::Evolution<Traits> experiment(pop, points);
  experiment.evolve(generations);

  Matrix R = experiment.get_best();
  check_equality(Traits::compute_determinant(R), 1);
}

void test_long_tetrahedron(const std::string fname)
{
  std::ifstream input(fname);
  CGAL::Surface_mesh<Point_3> mesh;
  if(!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file." << std::endl;
    std::exit(1);
  }

  // points in a matrix
  Matrix points;
  CGAL::Optimal_bounding_box::sm_to_matrix(mesh, points);

  std::size_t max_generations = 10;
  CGAL::Optimal_bounding_box::internal::Population<Traits> pop(50);
  CGAL::Optimal_bounding_box::internal::Evolution<Traits> experiment(pop, points);
  experiment.evolve(max_generations);

  Matrix R = experiment.get_best();
  check_equality(Traits::compute_determinant(R), 1);
  check_equality(R(0,0), -1);
  check_equality(R(0,1), 0);
  check_equality(R(0,2), 0);
  check_equality(R(1,0), 0);
  check_equality(R(1,1), -0.707107);
  check_equality(R(1,2), 0.707106) || assert_doubles(R(1,2), -0.707106);
  check_equality(R(2,0), 0);
  check_equality(R(2,1), 0.707106) || assert_doubles(R(1,2), -0.707106);
  check_equality(R(2,2), 0.707107);
}

void test_compute_obb_evolution(const std::string fname)
{
  std::ifstream input(fname);
  typedef CGAL::Surface_mesh<Point_3> SMesh;
  SMesh mesh;
  if(!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file." << std::endl;
    std::exit(1);
  }

  Traits traits;
  std::array<Point_3, 8> obb_points;
  CGAL::oriented_bounding_box(sm, obb_points, CGAL::parameters::use_convex_hull(true)
                                                               .geom_traits(traits));

  FT vol = CGAL::Optimal_bounding_box::calculate_volume(obb_points);
  check_equality(vol, 0.883371);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  CGAL::Surface_mesh<Point_3> result_mesh;
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
  CGAL::Surface_mesh<Point_3> mesh;
  if(!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file." << std::endl;
    std::exit(1);
  }

  CGAL::Surface_mesh<Point_3> obbmesh;
  CGAL::oriented_bounding_box(mesh, obbmesh);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_TEST
  std::ofstream out("/tmp/result_elephant.off");
  out << obbmesh;
  out.close();
#endif
}

void test_function_defaults_traits(const std::string fname1)
{
  std::ifstream input1(fname1);
  CGAL::Surface_mesh<Point_3> mesh1;
  if(!input1 || !(input1 >> mesh1) || mesh1.is_empty())
  {
    std::cerr << fname1 << " is not a valid off file." << std::endl;
    std::exit(1);
  }

  std::array<Point_3, 8> obb_points;
  CGAL::oriented_bounding_box(sm_points, obb_points, CGAL::parameters::use_convex_hull(true));

  const FT vol = CGAL::Optimal_bounding_box::calculate_volume(obb_points);
  check_equality(vol, 0.883371);
}

int main()
{
  test_genetic_algorithm();

  test_random_unit_tetra();
  test_reference_tetrahedron("data/reference_tetrahedron.off");
  test_long_tetrahedron("data/long_tetrahedron.off");
  test_compute_obb_evolution("data/random_unit_tetra.off");
  test_compute_obb_mesh("data/elephant.off");
  test_function_defaults_traits("data/random_unit_tetra.off");

  return 0;
}
