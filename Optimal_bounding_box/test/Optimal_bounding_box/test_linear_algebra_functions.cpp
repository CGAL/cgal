#include <CGAL/Optimal_bounding_box/linear_algebra.h>
#include <CGAL/Optimal_bounding_box/fitness_function.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>


typedef Eigen::Matrix3d Matrix3d;

bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}

void test_qr_factorization()
{
  Matrix3d A(3, 3);
  A << 0.3011944, 0.9932761, 0.5483701,
       0.5149142, 0.5973263, 0.5162336,
       0.0039213, 0.0202949, 0.9240308;

  Matrix3d Q;
  CGAL::Optimal_bounding_box::qr_factorization(A, Q);

  double epsilon = 1e-6;
  CGAL_assertion(assert_doubles(Q(0,0), -0.504895, epsilon));
  CGAL_assertion(assert_doubles(Q(0,1), 0.862834, epsilon));
  CGAL_assertion(assert_doubles(Q(0,2), -0.024447, epsilon));
  CGAL_assertion(assert_doubles(Q(1,0), -0.863156, epsilon));
  CGAL_assertion(assert_doubles(Q(1,1), -0.504894, epsilon));
  CGAL_assertion(assert_doubles(Q(1,2), 0.006687, epsilon));
  CGAL_assertion(assert_doubles(Q(2,0), -0.006573, epsilon));
  CGAL_assertion(assert_doubles(Q(2,1), 0.024478, epsilon));
  CGAL_assertion(assert_doubles(Q(2,2), 0.999679, epsilon));
}

void test_fitness_function()
{
  Eigen::Matrix<double, 4, 3> data_points(4, 3);
  data_points << 0.866802, 0.740808, 0.895304,
                 0.912651, 0.761565, 0.160330,
                 0.093661, 0.892578, 0.737412,
                 0.166461, 0.149912, 0.364944;

  Matrix3d rotation(3, 3);
  rotation << -0.809204, 0.124296, 0.574230,
             -0.574694, 0.035719, -0.817589,
             -0.122134, -0.991602, 0.042528;

  double fitness = CGAL::Optimal_bounding_box::compute_fitness(rotation, data_points);
  CGAL_assertion(assert_doubles(fitness, 0.58606, 1e-5));

}

void test_simplex_operations()
{
  Matrix3d Sc(3, 3);
  Sc << -0.809204, 0.124296, 0.574230,
        -0.574694, 0.035719, -0.817589,
        -0.122134, -0.991602, 0.042528;

  Matrix3d S_worst(3, 3);
  S_worst << -0.45070,  -0.32769,  -0.83035,
             -0.13619,  -0.89406,   0.42675,
             -0.88222,   0.30543,   0.35833;

  Matrix3d Sr = CGAL::Optimal_bounding_box::reflection(Sc, S_worst);
  double epsilon = 1e-5;
  CGAL_assertion(assert_doubles(Sr(0,0), -0.13359, epsilon));
  CGAL_assertion(assert_doubles(Sr(0,1), -0.95986, epsilon));
  CGAL_assertion(assert_doubles(Sr(0,2), -0.24664, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,0), -0.60307, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,1), -0.11875, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,2), 0.78880, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,0), -0.78642, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,1), 0.25411, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,2), -0.56300, epsilon));

  Matrix3d Se = CGAL::Optimal_bounding_box::expansion(Sc, S_worst, Sr);
  CGAL_assertion(assert_doubles(Se(0,0), -0.87991, epsilon));
  CGAL_assertion(assert_doubles(Se(0,1), 0.36105, epsilon));
  CGAL_assertion(assert_doubles(Se(0,2), -0.30888, epsilon));
  CGAL_assertion(assert_doubles(Se(1,0), -0.11816, epsilon));
  CGAL_assertion(assert_doubles(Se(1,1), -0.79593, epsilon));
  CGAL_assertion(assert_doubles(Se(1,2), -0.59375, epsilon));
  CGAL_assertion(assert_doubles(Se(2,0), -0.460215, epsilon));
  CGAL_assertion(assert_doubles(Se(2,1), -0.48595, epsilon));
  CGAL_assertion(assert_doubles(Se(2,2), 0.74300, epsilon));

  Matrix3d S_a(3, 3);
  S_a << -0.277970,  0.953559,  0.116010,
         -0.567497,  -0.065576,   -0.820760,
         -0.775035,   -0.293982,   0.559370;

  Matrix3d S_b(3, 3);
  S_b << -0.419979,   0.301765,  -0.855894,
         -0.653011,  -0.755415,   0.054087,
         -0.630234,   0.581624,   0.514314;

  Matrix3d S_c = CGAL::Optimal_bounding_box::mean(S_a, S_b);
  CGAL_assertion(assert_doubles(S_c(0,0), -0.35111, epsilon));
  CGAL_assertion(assert_doubles(S_c(0,1), 0.79308, epsilon));
  CGAL_assertion(assert_doubles(S_c(0,2), -0.49774, epsilon));
  CGAL_assertion(assert_doubles(S_c(1,0), -0.61398, epsilon));
  CGAL_assertion(assert_doubles(S_c(1,1), -0.59635, epsilon));
  CGAL_assertion(assert_doubles(S_c(1,2), -0.51710, epsilon));
  CGAL_assertion(assert_doubles(S_c(2,0), -0.70693, epsilon));
  CGAL_assertion(assert_doubles(S_c(2,1), 0.12405, epsilon));
  CGAL_assertion(assert_doubles(S_c(2,2), 0.69632, epsilon));
}

void test_centroid()
{
  Matrix3d S_a(3, 3);
  S_a <<  -0.588443,   0.807140,  -0.047542,
          -0.786228,  -0.584933,  -0.199246,
          -0.188629,  -0.079867,   0.978795;

  Matrix3d S_b(3, 3);
  S_b <<  -0.2192721,   0.2792986,  -0.9348326,
          -0.7772152,  -0.6292092,  -0.0056861,
          -0.5897934,   0.7253193,   0.3550431;

  Matrix3d S_c(3, 3);
  S_c <<   -0.32657,  -0.60013,  -0.73020,
           -0.20022,  -0.71110,   0.67398,
           -0.92372,   0.36630,   0.11207;

  Matrix3d S_centroid = CGAL::Optimal_bounding_box::centroid(S_a, S_b, S_c);
  double epsilon = 1e-5;
  CGAL_assertion(assert_doubles(S_centroid(0,0), -0.419979, epsilon));
  CGAL_assertion(assert_doubles(S_centroid(0,1), 0.301765, epsilon));
  CGAL_assertion(assert_doubles(S_centroid(0,2), -0.855894, epsilon));
  CGAL_assertion(assert_doubles(S_centroid(1,0), -0.653011, epsilon));
  CGAL_assertion(assert_doubles(S_centroid(1,1), -0.755415, epsilon));
  CGAL_assertion(assert_doubles(S_centroid(1,2), 0.054087, epsilon));
  CGAL_assertion(assert_doubles(S_centroid(2,0), -0.630234, epsilon));
  CGAL_assertion(assert_doubles(S_centroid(2,1), 0.581624, epsilon));
  CGAL_assertion(assert_doubles(S_centroid(2,2), 0.514314, epsilon));
}

int main()
{
  test_qr_factorization();
  test_fitness_function();
  test_simplex_operations();
  test_centroid();

  return 0;
}
