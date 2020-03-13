#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/optimal_bounding_box.h>

#include <CGAL/assertions.h>

#include <array>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point_3;

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

void test_fitness_function()
{
  std::array<Point_3, 4> points;
  points[0] = Point_3(0.866802, 0.740808, 0.895304);
  points[1] = Point_3(0.912651, 0.761565, 0.160330);
  points[2] = Point_3(0.093661, 0.892578, 0.737412);
  points[3] = Point_3(0.166461, 0.149912, 0.364944);

  Matrix rotation;
  rotation.set(0, 0, -0.809204);
  rotation.set(0, 1, 0.124296);
  rotation.set(0, 2, 0.574230);
  rotation.set(1, 0, -0.574694);
  rotation.set(1, 1, 0.035719);
  rotation.set(1, 2, -0.817589);
  rotation.set(2, 0, -0.122134);
  rotation.set(2, 1, -0.991602);
  rotation.set(2, 2, 0.042528);

  const double fitness = CGAL::Optimal_bounding_box::internal::compute_fitness<Traits>(rotation, points);
  check_equality(fitness, 0.58606);
}

void test_eigen_matrix_interface()
{
  Matrix A;
  A.set(0, 0, 0.1);
  A.set(0, 1, 0.2);
  A.set(0, 2, 0.3);
  A.set(1, 0, 0.4);
  A.set(1, 1, 0.5);
  A.set(1, 2, 0.6);
  A.set(2, 0, 0.7);
  A.set(2, 1, 0.8);
  A.set(2, 2, 0.9);

  Matrix B;
  B = CGAL::Optimal_bounding_box::internal::transpose(A);

  Matrix S;
  S = 0.5 * A;

  Matrix C;
  C.set(0, 0, 0.3011944);
  C.set(0, 1, 0.9932761);
  C.set(0, 2, 0.5483701);
  C.set(1, 0, 0.5149142);
  C.set(1, 1, 0.5973263);
  C.set(1, 2, 0.5162336);
  C.set(2, 0, 0.0039213);
  C.set(2, 1, 0.0202949);
  C.set(2, 2, 0.9240308);

  Matrix Q = Traits::get_Q(C);

  check_equality(Q(0,0), -0.504895);
  check_equality(Q(0,1), 0.862834);
  check_equality(Q(0,2), -0.024447);
  check_equality(Q(1,0), -0.863156);
  check_equality(Q(1,1), -0.504894);
  check_equality(Q(1,2), 0.006687);
  check_equality(Q(2,0), -0.006573);
  check_equality(Q(2,1), 0.024478);
  check_equality(Q(2,2), 0.999679);

  Matrix D;
  D.set(0, 0, -0.809204);
  D.set(0, 1, 0.124296);
  D.set(0, 2, 0.574230);
  D.set(1, 0, -0.574694);
  D.set(1, 1, 0.035719);
  D.set(1, 2, -0.817589);
  D.set(2, 0, -0.122134);
  D.set(2, 1, -0.991602);
  D.set(2, 2, 0.042528);

  Matrix E;
  E.set(0, 0, -0.45070);
  E.set(0, 1, -0.32769);
  E.set(0, 2, -0.83035);
  E.set(1, 0, -0.13619);
  E.set(1, 1, -0.89406);
  E.set(1, 2, 0.42675);
  E.set(2, 0, -0.88222);
  E.set(2, 1, 0.30543);
  E.set(2, 2, 0.35833);

  Matrix Sr = CGAL::Optimal_bounding_box::internal::reflection(D, E);
  check_equality(Sr(0,0), -0.13359);
  check_equality(Sr(0,1), -0.95986);
  check_equality(Sr(0,2), -0.24664);
  check_equality(Sr(1,0), -0.60307);
  check_equality(Sr(1,1), -0.11875);
  check_equality(Sr(1,2), 0.78880);
  check_equality(Sr(2,0), -0.78642);
  check_equality(Sr(2,1), 0.25411);
  check_equality(Sr(2,2), -0.56300);
}

int main(int, char**)
{
  test_fitness_function();
  test_eigen_matrix_interface();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
