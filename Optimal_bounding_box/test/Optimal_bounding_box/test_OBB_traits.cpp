#include <CGAL/Eigen_linear_algebra_traits.h>
#include <CGAL/Optimal_bounding_box/nelder_mead_functions.h>
#include <CGAL/Optimal_bounding_box/fitness_function.h>

#include <CGAL/assertions.h>

bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
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
  test_fitness_function();
  test_eigen_matrix_interface();

  return 0;
}
