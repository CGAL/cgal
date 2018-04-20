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
  CGAL_assertion(assert_doubles(fitness, 0.58606, 1e-6));

}






int main()
{

  test_qr_factorization();
  test_fitness_function();

  return 0;
}
