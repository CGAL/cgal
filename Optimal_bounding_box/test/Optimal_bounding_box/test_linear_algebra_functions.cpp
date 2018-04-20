#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Optimal_bounding_box/linear_algebra.h>
#include <CGAL/Optimal_bounding_box/fitness_function.h>

#include <iostream>
#include <fstream>

#include <Eigen/Dense>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;

typedef Eigen::MatrixXf MatrixXf;


bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}


void test_qr_factorization()
{

  MatrixXf A(3, 3);
  A << 0.3011944, 0.9932761, 0.5483701,
       0.5149142, 0.5973263, 0.5162336,
       0.0039213, 0.0202949, 0.9240308;

  MatrixXf Q;
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
  MatrixXf data_points(4, 3);
  data_points << 0.866802, 0.740808, 0.895304,
                 0.912651, 0.761565, 0.160330,
                 0.093661, 0.892578, 0.737412,
                 0.166461, 0.149912, 0.364944;

  MatrixXf rotation(3, 3);
  rotation << -0.809204, 0.124296, 0.574230,
             -0.574694, 0.035719, -0.817589,
             -0.122134, -0.991602, 0.042528;

  double fitness = CGAL::Optimal_bounding_box::compute_fitness(rotation, data_points);
  CGAL_assertion(assert_doubles(fitness, 0.58606, 1e-5));

}


void test_simplex_operations()
{
  MatrixXf Sc(3, 3);
  Sc << -0.809204, 0.124296, 0.574230,
        -0.574694, 0.035719, -0.817589,
        -0.122134, -0.991602, 0.042528;

  MatrixXf S_worst(3, 3);
  S_worst << -0.45070,  -0.32769,  -0.83035,
             -0.13619,  -0.89406,   0.42675,
             -0.88222,   0.30543,   0.35833;

  MatrixXf Sr = CGAL::Optimal_bounding_box::reflection(Sc, S_worst);
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

  MatrixXf Se = CGAL::Optimal_bounding_box::expansion(Sc, S_worst, Sr);
  CGAL_assertion(assert_doubles(Se(0,0), -0.87991, epsilon));
  CGAL_assertion(assert_doubles(Se(0,1), 0.36105, epsilon));
  CGAL_assertion(assert_doubles(Se(0,2), -0.30888, epsilon));
  CGAL_assertion(assert_doubles(Se(1,0), -0.11816, epsilon));
  CGAL_assertion(assert_doubles(Se(1,1), -0.79593, epsilon));
  CGAL_assertion(assert_doubles(Se(1,2), -0.59375, epsilon));
  CGAL_assertion(assert_doubles(Se(2,0), -0.460215, epsilon));
  CGAL_assertion(assert_doubles(Se(2,1), -0.48595, epsilon));
  CGAL_assertion(assert_doubles(Se(2,2), 0.74300, epsilon));


  MatrixXf S_a(3, 3);
  S_a << -0.277970,  0.953559,  0.116010,
         -0.567497,  -0.065576,   -0.820760,
         -0.775035,   -0.293982,   0.559370;

  MatrixXf S_b(3, 3);
  S_b << -0.419979,   0.301765,  -0.855894,
         -0.653011,  -0.755415,   0.054087,
         -0.630234,   0.581624,   0.514314;

  MatrixXf S_c = CGAL::Optimal_bounding_box::constraction(S_a, S_b);
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





int main()
{

  test_qr_factorization();
  test_fitness_function();
  test_simplex_operations();

  return 0;
}
