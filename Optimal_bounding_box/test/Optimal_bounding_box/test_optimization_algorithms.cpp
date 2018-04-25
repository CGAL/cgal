#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Optimal_bounding_box/optimization_algorithms.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/obb.h>

#include <iostream>
#include <fstream>

#include <Eigen/Dense>


typedef Eigen::MatrixXf MatrixXf;

bool assert_doubles(double d1, double d2, double epsilon)
{
  return (d1 < d2 + epsilon && d1 > d2 - epsilon) ? true : false;
}


void test_population()
{

  CGAL::Optimal_bounding_box::Population<MatrixXf> pop(5);

  //pop.show_population();

  CGAL_assertion(pop.size() == 5);
}

void test_nelder_mead()
{

  MatrixXf data_points(4, 3);
  data_points << 0.866802, 0.740808, 0.895304,
                 0.912651, 0.761565, 0.160330,
                 0.093661, 0.892578, 0.737412,
                 0.166461, 0.149912, 0.364944;

  // one simplex
  std::vector<MatrixXf> simplex(4);

  MatrixXf v0(3, 3);
  v0 <<  -0.2192721,   0.2792986,  -0.9348326,
          -0.7772152,  -0.6292092,  -0.0056861,
          -0.5897934,   0.7253193,   0.3550431;

  MatrixXf v1(3, 3);
  v1 <<  -0.588443,   0.807140,  -0.047542,
          -0.786228,  -0.584933,  -0.199246,
          -0.188629,  -0.079867,   0.978795;

  MatrixXf v2(3, 3);
  v2 << -0.277970,  0.953559,  0.116010,
         -0.567497,  -0.065576,   -0.820760,
         -0.775035,   -0.293982,   0.559370;

  MatrixXf v3(3, 3);
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

  MatrixXf v0_new = simplex[0];
  CGAL_assertion(assert_doubles(v0_new(0,0), -0.288975, epsilon));
  CGAL_assertion(assert_doubles(v0_new(0,1), 0.7897657, epsilon));
  CGAL_assertion(assert_doubles(v0_new(0,2), -0.541076, epsilon));
  CGAL_assertion(assert_doubles(v0_new(1,0), -0.9407046, epsilon));
  CGAL_assertion(assert_doubles(v0_new(1,1), -0.3391466, epsilon));
  CGAL_assertion(assert_doubles(v0_new(1,2), 0.0073817, epsilon));
  CGAL_assertion(assert_doubles(v0_new(2,0), -0.1776743, epsilon));
  CGAL_assertion(assert_doubles(v0_new(2,1), 0.5111260, epsilon));
  CGAL_assertion(assert_doubles(v0_new(2,2), 0.84094, epsilon));

  MatrixXf v1_new = simplex[1];
  CGAL_assertion(assert_doubles(v1_new(0,0), -0.458749, epsilon));
  CGAL_assertion(assert_doubles(v1_new(0,1), 0.823283, epsilon));
  CGAL_assertion(assert_doubles(v1_new(0,2), -0.334296, epsilon));
  CGAL_assertion(assert_doubles(v1_new(1,0), -0.885235, epsilon));
  CGAL_assertion(assert_doubles(v1_new(1,1), -0.455997, epsilon));
  CGAL_assertion(assert_doubles(v1_new(1,2), 0.091794, epsilon));
  CGAL_assertion(assert_doubles(v1_new(2,0), -0.076866, epsilon));
  CGAL_assertion(assert_doubles(v1_new(2,1), 0.338040, epsilon));
  CGAL_assertion(assert_doubles(v1_new(2,2), 0.937987, epsilon));

  MatrixXf v2_new = simplex[2];
  CGAL_assertion(assert_doubles(v2_new(0,0), -0.346582, epsilon));
  CGAL_assertion(assert_doubles(v2_new(0,1), 0.878534, epsilon));
  CGAL_assertion(assert_doubles(v2_new(0,2), -0.328724, epsilon));
  CGAL_assertion(assert_doubles(v2_new(1,0), -0.936885, epsilon));
  CGAL_assertion(assert_doubles(v2_new(1,1), -0.341445, epsilon));
  CGAL_assertion(assert_doubles(v2_new(1,2), 0.075251, epsilon));
  CGAL_assertion(assert_doubles(v2_new(2,0), -0.046131, epsilon));
  CGAL_assertion(assert_doubles(v2_new(2,1), 0.334057, epsilon));
  CGAL_assertion(assert_doubles(v2_new(2,2), 0.941423, epsilon));

  MatrixXf v3_new = simplex[3];
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

  MatrixXf data_points(4, 3);
  data_points << 0.866802, 0.740808, 0.895304,
                 0.912651, 0.761565, 0.160330,
                 0.093661, 0.892578, 0.737412,
                 0.166461, 0.149912, 0.364944;

  CGAL::Optimal_bounding_box::Population<MatrixXf> pop(5);
  CGAL::Optimal_bounding_box::genetic_algorithm(pop, data_points);
  CGAL_assertion(pop.size() == 5);
}

void test_visualization()
{
  MatrixXf data_points(4, 3);
  data_points << 0.866802, 0.740808, 0.895304,
                 0.912651, 0.761565, 0.160330,
                 0.093661, 0.892578, 0.737412,
                 0.166461, 0.149912, 0.364944;


  CGAL::Optimal_bounding_box::visualize_obb(data_points);


}


int main()
{
  test_population();
  test_nelder_mead();
  test_genetic_algorithm();


  test_visualization();







  return 0;
}
