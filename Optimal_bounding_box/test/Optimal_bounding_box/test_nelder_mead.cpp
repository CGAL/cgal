#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/optimal_bounding_box.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Point_3                                                  Point_3;

typedef CGAL::Oriented_bounding_box_traits<K>                       Traits;
typedef Traits::Matrix                                              Matrix;

bool are_equal(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}

void test_simplex_operations()
{
  const double epsilon = 1e-5;

  Matrix Sc(3, 3);
  Sc.set(0, 0, -0.809204); Sc.set(0, 1, 0.124296); Sc.set(0, 2, 0.574230);
  Sc.set(1, 0, -0.574694); Sc.set(1, 1, 0.035719); Sc.set(1, 2, -0.817589);
  Sc.set(2, 0, -0.122134); Sc.set(2, 1, -0.991602); Sc.set(2, 2, 0.042528);

  Matrix S_worst(3, 3);
  S_worst.set(0, 0, -0.45070); S_worst.set(0, 1, -0.32769); S_worst.set(0, 2, -0.83035);
  S_worst.set(1, 0, -0.13619); S_worst.set(1, 1, -0.89406); S_worst.set(1, 2, 0.42675);
  S_worst.set(2, 0, -0.88222); S_worst.set(2, 1, 0.30543); S_worst.set(2, 2, 0.35833);

  Matrix Sr = CGAL::Optimal_bounding_box::reflection<Linear_algebra_traits>(Sc, S_worst);
  assert(are_equal(Sr(0,0), -0.13359, epsilon));
  assert(are_equal(Sr(0,1), -0.95986, epsilon));
  assert(are_equal(Sr(0,2), -0.24664, epsilon));
  assert(are_equal(Sr(1,0), -0.60307, epsilon));
  assert(are_equal(Sr(1,1), -0.11875, epsilon));
  assert(are_equal(Sr(1,2), 0.78880, epsilon));
  assert(are_equal(Sr(2,0), -0.78642, epsilon));
  assert(are_equal(Sr(2,1), 0.25411, epsilon));
  assert(are_equal(Sr(2,2), -0.56300, epsilon));

  Matrix Se = CGAL::Optimal_bounding_box::expansion<Linear_algebra_traits>(Sc, S_worst, Sr);
  assert(are_equal(Se(0,0), -0.87991, epsilon));
  assert(are_equal(Se(0,1), 0.36105, epsilon));
  assert(are_equal(Se(0,2), -0.30888, epsilon));
  assert(are_equal(Se(1,0), -0.11816, epsilon));
  assert(are_equal(Se(1,1), -0.79593, epsilon));
  assert(are_equal(Se(1,2), -0.59375, epsilon));
  assert(are_equal(Se(2,0), -0.460215, epsilon));
  assert(are_equal(Se(2,1), -0.48595, epsilon));
  assert(are_equal(Se(2,2), 0.74300, epsilon));

  Matrix S_a(3, 3);
  S_a.set(0, 0, -0.277970); S_a.set(0, 1, 0.953559); S_a.set(0, 2, 0.116010);
  S_a.set(1, 0, -0.567497); S_a.set(1, 1, -0.065576); S_a.set(1, 2, -0.820760);
  S_a.set(2, 0, -0.775035); S_a.set(2, 1, -0.293982); S_a.set(2, 2, 0.559370);

  Matrix S_b(3, 3);
  S_b.set(0, 0, -0.419979); S_b.set(0, 1, 0.301765); S_b.set(0, 2, -0.8558940);
  S_b.set(1, 0, -0.653011); S_b.set(1, 1, -0.755415); S_b.set(1, 2, 0.054087);
  S_b.set(2, 0, -0.630234); S_b.set(2, 1, 0.581624); S_b.set(2, 2, 0.514314);

  Matrix S_c = CGAL::Optimal_bounding_box::mean<Linear_algebra_traits>(S_a, S_b);
  assert(are_equal(S_c(0,0), -0.35111, epsilon));
  assert(are_equal(S_c(0,1), 0.79308, epsilon));
  assert(are_equal(S_c(0,2), -0.49774, epsilon));
  assert(are_equal(S_c(1,0), -0.61398, epsilon));
  assert(are_equal(S_c(1,1), -0.59635, epsilon));
  assert(are_equal(S_c(1,2), -0.51710, epsilon));
  assert(are_equal(S_c(2,0), -0.70693, epsilon));
  assert(are_equal(S_c(2,1), 0.12405, epsilon));
  assert(are_equal(S_c(2,2), 0.69632, epsilon));
}

void test_centroid()
{
  const double epsilon = 1e-5;

  Matrix S_a;
  S_a.set(0, 0, -0.588443); S_a.set(0, 1, 0.807140); S_a.set(0, 2, -0.047542);
  S_a.set(1, 0, -0.786228); S_a.set(1, 1, -0.584933); S_a.set(1, 2, -0.199246);
  S_a.set(2, 0, -0.188629); S_a.set(2, 1, -0.079867); S_a.set(2, 2, 0.978795);

  Matrix S_b(3, 3);
  S_b.set(0, 0, -0.2192721); S_b.set(0, 1, 0.2792986); S_b.set(0, 2, -0.9348326);
  S_b.set(1, 0, -0.7772152); S_b.set(1, 1, -0.6292092); S_b.set(1, 2, -0.005686);
  S_b.set(2, 0, -0.5897934); S_b.set(2, 1, 0.7253193); S_b.set(2, 2, 0.3550431);

  Matrix S_c(3, 3);
  S_c.set(0, 0, -0.32657); S_c.set(0, 1, -0.60013); S_c.set(0, 2, -0.730206);
  S_c.set(1, 0, -0.20022); S_c.set(1, 1, -0.71110); S_c.set(1, 2, 0.67398);
  S_c.set(2, 0, -0.92372); S_c.set(2, 1, 0.36630); S_c.set(2, 2, 0.11207);

  Matrix S_centroid = CGAL::Optimal_bounding_box::nm_centroid<Linear_algebra_traits>(S_a, S_b, S_c);
  assert(are_equal(S_centroid(0,0), -0.419979, epsilon));
  assert(are_equal(S_centroid(0,1), 0.301765, epsilon));
  assert(are_equal(S_centroid(0,2), -0.855894, epsilon));
  assert(are_equal(S_centroid(1,0), -0.653011, epsilon));
  assert(are_equal(S_centroid(1,1), -0.755415, epsilon));
  assert(are_equal(S_centroid(1,2), 0.054087, epsilon));
  assert(are_equal(S_centroid(2,0), -0.630234, epsilon));
  assert(are_equal(S_centroid(2,1), 0.581624, epsilon));
  assert(are_equal(S_centroid(2,2), 0.514314, epsilon));
}

void test_nelder_mead()
{
  std::array<Point_3, 4> points;
  points[0] = Point_3(0.866802, 0.740808, 0.895304);
  points[1] = Point_3(0.912651, 0.761565, 0.160330);
  points[2] = Point_3(0.093661, 0.892578, 0.737412);
  points[3] = Point_3(0.166461, 0.149912, 0.364944);

  // one simplex
  std::array<Matrix, 4> simplex;

  Matrix v0, v1, v2, v3;
  v0(0,0) = -0.2192721; v0(0,1) = 0.2792986; v0(0,2) = -0.9348326;
  v0(1,0) = -0.7772152; v0(1,1) = -0.6292092; v0(1,2) = -0.0056861;
  v0(2,0) = -0.5897934; v0(2,1) = 0.7253193; v0(2,2) = 0.3550431;

  v1(0,0) = -0.588443; v1(0,1) = 0.807140; v1(0,2) = -0.047542;
  v1(1,0) = -0.786228; v1(1,1) = -0.584933; v1(1,2) = -0.199246;
  v1(2,0) = -0.188629; v1(2,1) = -0.079867; v1(2,2) = 0.978795;

  v2(0,0) = -0.277970; v2(0,1) = 0.953559; v2(0,2) = 0.116010;
  v2(1,0) = -0.567497; v2(1,1) = -0.065576; v2(1,2) = -0.820760;
  v2(2,0) = -0.775035; v2(2,1) = -0.293982; v2(2,2) = 0.559370;

  v3(0,0) = -0.32657; v3(0,1) = -0.60013; v3(0,2) = -0.73020;
  v3(1,0) = -0.20022; v3(1,1) = -0.71110; v3(1,2) = 0.67398;
  v3(2,0) = -0.92372; v3(2,1) = 0.36630; v3(2,2) = 0.11207;

  simplex[0] = v0;
  simplex[1] = v1;
  simplex[2] = v2;
  simplex[3] = v3;

  std::size_t nm_iterations = 19;
  CGAL::Optimal_bounding_box::internal::nelder_mead<Linear_algebra_traits>(simplex, data_points, nm_iterations);

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

int main()
{
  test_simplex_operations();
  test_centroid();
  test_nelder_mead();

  return 0;
}
