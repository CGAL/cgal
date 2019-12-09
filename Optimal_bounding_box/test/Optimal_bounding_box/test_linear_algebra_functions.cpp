#include <CGAL/Eigen_linear_algebra_traits.h>
#include <CGAL/Optimal_bounding_box/nelder_mead_functions.h>
#include <CGAL/Optimal_bounding_box/fitness_function.h>

#include <CGAL/assertions.h>

bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}

void test_get_Q()
{
  typedef CGAL::Eigen_dense_matrix<double, 3, 3> Mat;
  Mat A(3, 3);
  A.set_coef(0, 0, 0.3011944);
  A.set_coef(0, 1, 0.9932761);
  A.set_coef(0, 2, 0.5483701);
  A.set_coef(1, 0, 0.5149142);
  A.set_coef(1, 1, 0.5973263);
  A.set_coef(1, 2, 0.5162336);
  A.set_coef(2, 0, 0.0039213);
  A.set_coef(2, 1, 0.0202949);
  A.set_coef(2, 2, 0.9240308);

  CGAL_assertion_code(Mat Q = CGAL::Eigen_linear_algebra_traits::get_Q(A));
  CGAL_assertion_code(double epsilon = 1e-6);
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
  typedef typename CGAL::Eigen_linear_algebra_traits::MatrixXd MatrixXd;
  MatrixXd data_points(4, 3);

  data_points.set_coef(0, 0, 0.866802);
  data_points.set_coef(0, 1, 0.740808);
  data_points.set_coef(0, 2, 0.895304);

  data_points.set_coef(1, 0, 0.912651);
  data_points.set_coef(1, 1, 0.761565);
  data_points.set_coef(1, 2, 0.160330);

  data_points.set_coef(2, 0, 0.093661);
  data_points.set_coef(2, 1, 0.892578);
  data_points.set_coef(2, 2, 0.737412);

  data_points.set_coef(3, 0, 0.166461);
  data_points.set_coef(3, 1, 0.149912);
  data_points.set_coef(3, 2, 0.364944);

  typedef typename CGAL::Eigen_linear_algebra_traits::Matrix3d Matrix3d;
  Matrix3d rotation;
  rotation.set_coef(0, 0, -0.809204);
  rotation.set_coef(0, 1, 0.124296);
  rotation.set_coef(0, 2, 0.574230);
  rotation.set_coef(1, 0, -0.574694);
  rotation.set_coef(1, 1, 0.035719);
  rotation.set_coef(1, 2, -0.817589);
  rotation.set_coef(2, 0, -0.122134);
  rotation.set_coef(2, 1, -0.991602);
  rotation.set_coef(2, 2, 0.042528);

  CGAL_assertion_code(typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits);
  CGAL_assertion_code(double fitness = CGAL::Optimal_bounding_box::
      compute_fitness<Linear_algebra_traits> (rotation, data_points));
  CGAL_assertion(assert_doubles(fitness, 0.58606, 1e-5));
}

void test_simplex_operations()
{
  CGAL_assertion_code(typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits);
  typedef CGAL::Eigen_dense_matrix<double, 3, 3> Matrix;

  Matrix Sc(3, 3);
  Sc.set_coef(0, 0, -0.809204);
  Sc.set_coef(0, 1, 0.124296);
  Sc.set_coef(0, 2, 0.574230);
  Sc.set_coef(1, 0, -0.574694);
  Sc.set_coef(1, 1, 0.035719);
  Sc.set_coef(1, 2, -0.817589);
  Sc.set_coef(2, 0, -0.122134);
  Sc.set_coef(2, 1, -0.991602);
  Sc.set_coef(2, 2, 0.042528);

  Matrix S_worst(3, 3);
  S_worst.set_coef(0, 0, -0.45070);
  S_worst.set_coef(0, 1, -0.32769);
  S_worst.set_coef(0, 2, -0.83035);
  S_worst.set_coef(1, 0, -0.13619);
  S_worst.set_coef(1, 1, -0.89406);
  S_worst.set_coef(1, 2, 0.42675);
  S_worst.set_coef(2, 0, -0.88222);
  S_worst.set_coef(2, 1, 0.30543);
  S_worst.set_coef(2, 2, 0.35833);

  CGAL_assertion_code(Matrix Sr = CGAL::Optimal_bounding_box::reflection<Linear_algebra_traits>(Sc, S_worst));
  CGAL_assertion_code(double epsilon = 1e-5);
  CGAL_assertion(assert_doubles(Sr(0,0), -0.13359, epsilon));
  CGAL_assertion(assert_doubles(Sr(0,1), -0.95986, epsilon));
  CGAL_assertion(assert_doubles(Sr(0,2), -0.24664, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,0), -0.60307, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,1), -0.11875, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,2), 0.78880, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,0), -0.78642, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,1), 0.25411, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,2), -0.56300, epsilon));

  CGAL_assertion_code(Matrix Se =
      CGAL::Optimal_bounding_box::expansion<Linear_algebra_traits>(Sc, S_worst, Sr));
  CGAL_assertion(assert_doubles(Se(0,0), -0.87991, epsilon));
  CGAL_assertion(assert_doubles(Se(0,1), 0.36105, epsilon));
  CGAL_assertion(assert_doubles(Se(0,2), -0.30888, epsilon));
  CGAL_assertion(assert_doubles(Se(1,0), -0.11816, epsilon));
  CGAL_assertion(assert_doubles(Se(1,1), -0.79593, epsilon));
  CGAL_assertion(assert_doubles(Se(1,2), -0.59375, epsilon));
  CGAL_assertion(assert_doubles(Se(2,0), -0.460215, epsilon));
  CGAL_assertion(assert_doubles(Se(2,1), -0.48595, epsilon));
  CGAL_assertion(assert_doubles(Se(2,2), 0.74300, epsilon));

  Matrix S_a(3, 3);
  S_a.set_coef(0, 0, -0.277970);
  S_a.set_coef(0, 1, 0.953559);
  S_a.set_coef(0, 2, 0.116010);
  S_a.set_coef(1, 0, -0.567497);
  S_a.set_coef(1, 1, -0.065576);
  S_a.set_coef(1, 2, -0.820760);
  S_a.set_coef(2, 0, -0.775035);
  S_a.set_coef(2, 1, -0.293982);
  S_a.set_coef(2, 2, 0.559370);

  Matrix S_b(3, 3);
  S_b.set_coef(0, 0, -0.419979);
  S_b.set_coef(0, 1, 0.301765);
  S_b.set_coef(0, 2, -0.8558940);
  S_b.set_coef(1, 0, -0.653011);
  S_b.set_coef(1, 1, -0.755415);
  S_b.set_coef(1, 2, 0.054087);
  S_b.set_coef(2, 0, -0.630234);
  S_b.set_coef(2, 1, 0.581624);
  S_b.set_coef(2, 2, 0.514314);

  CGAL_assertion_code(Matrix S_c =
      CGAL::Optimal_bounding_box::mean<Linear_algebra_traits>(S_a, S_b));
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
  CGAL_assertion_code(typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits);
  typedef CGAL::Eigen_dense_matrix<double, 3, 3> Matrix;

  Matrix S_a;
  S_a.set_coef(0, 0, -0.588443);
  S_a.set_coef(0, 1, 0.807140);
  S_a.set_coef(0, 2, -0.047542);
  S_a.set_coef(1, 0, -0.786228);
  S_a.set_coef(1, 1, -0.584933);
  S_a.set_coef(1, 2, -0.199246);
  S_a.set_coef(2, 0, -0.188629);
  S_a.set_coef(2, 1, -0.079867);
  S_a.set_coef(2, 2, 0.978795);

  Matrix S_b(3, 3);
  S_b.set_coef(0, 0, -0.2192721);
  S_b.set_coef(0, 1, 0.2792986);
  S_b.set_coef(0, 2, -0.9348326);
  S_b.set_coef(1, 0, -0.7772152);
  S_b.set_coef(1, 1, -0.6292092);
  S_b.set_coef(1, 2, -0.005686);
  S_b.set_coef(2, 0, -0.5897934);
  S_b.set_coef(2, 1, 0.7253193);
  S_b.set_coef(2, 2, 0.3550431);

  Matrix S_c(3, 3);
  S_c.set_coef(0, 0, -0.32657);
  S_c.set_coef(0, 1, -0.60013);
  S_c.set_coef(0, 2, -0.730206);
  S_c.set_coef(1, 0, -0.20022);
  S_c.set_coef(1, 1, -0.71110);
  S_c.set_coef(1, 2, 0.67398);
  S_c.set_coef(2, 0, -0.92372);
  S_c.set_coef(2, 1, 0.36630);
  S_c.set_coef(2, 2, 0.11207);

  CGAL_assertion_code(Matrix S_centroid =
      CGAL::Optimal_bounding_box::nm_centroid<Linear_algebra_traits>(S_a, S_b, S_c));
  CGAL_assertion_code(double epsilon = 1e-5);
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

void test_eigen_matrix_interface()
{
  CGAL_assertion_code(typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits);
  typedef CGAL::Eigen_dense_matrix<double, 3, 3> Matrix;

  Matrix A(3, 3);
  A.set_coef(0, 0, 0.1);
  A.set_coef(0, 1, 0.2);
  A.set_coef(0, 2, 0.3);
  A.set_coef(1, 0, 0.4);
  A.set_coef(1, 1, 0.5);
  A.set_coef(1, 2, 0.6);
  A.set_coef(2, 0, 0.7);
  A.set_coef(2, 1, 0.8);
  A.set_coef(2, 2, 0.9);
  CGAL_assertion_code(Matrix B);
  CGAL_assertion_code(B = CGAL::Eigen_linear_algebra_traits::transpose(A));
  CGAL_assertion_code(Matrix S);
  CGAL_assertion_code(S = 0.5 * A);
  Matrix C(3,3);
  C.set_coef(0, 0, 0.3011944);
  C.set_coef(0, 1, 0.9932761);
  C.set_coef(0, 2, 0.5483701);
  C.set_coef(1, 0, 0.5149142);
  C.set_coef(1, 1, 0.5973263);
  C.set_coef(1, 2, 0.5162336);
  C.set_coef(2, 0, 0.0039213);
  C.set_coef(2, 1, 0.0202949);
  C.set_coef(2, 2, 0.9240308);

  CGAL_assertion_code(Matrix Q = CGAL::Eigen_linear_algebra_traits::get_Q(C));
  CGAL_assertion_code(double epsilon = 1e-5);
  CGAL_assertion(assert_doubles(Q(0,0), -0.504895, epsilon));
  CGAL_assertion(assert_doubles(Q(0,1), 0.862834, epsilon));
  CGAL_assertion(assert_doubles(Q(0,2), -0.024447, epsilon));
  CGAL_assertion(assert_doubles(Q(1,0), -0.863156, epsilon));
  CGAL_assertion(assert_doubles(Q(1,1), -0.504894, epsilon));
  CGAL_assertion(assert_doubles(Q(1,2), 0.006687, epsilon));
  CGAL_assertion(assert_doubles(Q(2,0), -0.006573, epsilon));
  CGAL_assertion(assert_doubles(Q(2,1), 0.024478, epsilon));
  CGAL_assertion(assert_doubles(Q(2,2), 0.999679, epsilon));

  Matrix D(3,3);
  D.set_coef(0, 0, -0.809204);
  D.set_coef(0, 1, 0.124296);
  D.set_coef(0, 2, 0.574230);
  D.set_coef(1, 0, -0.574694);
  D.set_coef(1, 1, 0.035719);
  D.set_coef(1, 2, -0.817589);
  D.set_coef(2, 0, -0.122134);
  D.set_coef(2, 1, -0.991602);
  D.set_coef(2, 2, 0.042528);

  Matrix E(3,3);
  E.set_coef(0, 0, -0.45070);
  E.set_coef(0, 1, -0.32769);
  E.set_coef(0, 2, -0.83035);
  E.set_coef(1, 0, -0.13619);
  E.set_coef(1, 1, -0.89406);
  E.set_coef(1, 2, 0.42675);
  E.set_coef(2, 0, -0.88222);
  E.set_coef(2, 1, 0.30543);
  E.set_coef(2, 2, 0.35833);

  CGAL_assertion_code(Matrix Sr = CGAL::Optimal_bounding_box::reflection<Linear_algebra_traits>(D, E));
  CGAL_assertion(assert_doubles(Sr(0,0), -0.13359, epsilon));
  CGAL_assertion(assert_doubles(Sr(0,1), -0.95986, epsilon));
  CGAL_assertion(assert_doubles(Sr(0,2), -0.24664, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,0), -0.60307, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,1), -0.11875, epsilon));
  CGAL_assertion(assert_doubles(Sr(1,2), 0.78880, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,0), -0.78642, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,1), 0.25411, epsilon));
  CGAL_assertion(assert_doubles(Sr(2,2), -0.56300, epsilon));
}

int main()
{
  test_get_Q();
  test_fitness_function();
  test_simplex_operations();
  test_centroid();
  test_eigen_matrix_interface();

  return 0;
}
