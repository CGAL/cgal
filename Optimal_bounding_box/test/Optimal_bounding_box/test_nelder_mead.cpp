#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/optimal_bounding_box.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point_3;

typedef CGAL::Oriented_bounding_box_traits_3<K>                     Traits;
typedef Traits::Matrix                                              Matrix;
typedef CGAL::Optimal_bounding_box::internal::Population<Traits>::Vertex Vertex;

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

void test_simplex_operations(const Traits& traits)
{
  Matrix Sc;
  Sc.set(0, 0, -0.809204); Sc.set(0, 1, 0.124296); Sc.set(0, 2, 0.574230);
  Sc.set(1, 0, -0.574694); Sc.set(1, 1, 0.035719); Sc.set(1, 2, -0.817589);
  Sc.set(2, 0, -0.122134); Sc.set(2, 1, -0.991602); Sc.set(2, 2, 0.042528);

  Matrix S_worst;
  S_worst.set(0, 0, -0.45070); S_worst.set(0, 1, -0.32769); S_worst.set(0, 2, -0.83035);
  S_worst.set(1, 0, -0.13619); S_worst.set(1, 1, -0.89406); S_worst.set(1, 2, 0.42675);
  S_worst.set(2, 0, -0.88222); S_worst.set(2, 1, 0.30543); S_worst.set(2, 2, 0.35833);

  Matrix Sr = CGAL::Optimal_bounding_box::internal::reflection(Sc, S_worst);
  check_equality(Sr(0,0), -0.13359);
  check_equality(Sr(0,1), -0.95986);
  check_equality(Sr(0,2), -0.24664);
  check_equality(Sr(1,0), -0.60307);
  check_equality(Sr(1,1), -0.11875);
  check_equality(Sr(1,2), 0.78880);
  check_equality(Sr(2,0), -0.78642);
  check_equality(Sr(2,1), 0.25411);
  check_equality(Sr(2,2), -0.56300);

  Matrix Se = CGAL::Optimal_bounding_box::internal::expansion(Sc, S_worst, Sr);
  check_equality(Se(0,0), -0.87991);
  check_equality(Se(0,1), 0.36105);
  check_equality(Se(0,2), -0.30888);
  check_equality(Se(1,0), -0.11816);
  check_equality(Se(1,1), -0.79593);
  check_equality(Se(1,2), -0.59375);
  check_equality(Se(2,0), -0.460215);
  check_equality(Se(2,1), -0.48595);
  check_equality(Se(2,2), 0.74300);

  Matrix S_a;
  S_a.set(0, 0, -0.277970); S_a.set(0, 1, 0.953559); S_a.set(0, 2, 0.116010);
  S_a.set(1, 0, -0.567497); S_a.set(1, 1, -0.065576); S_a.set(1, 2, -0.820760);
  S_a.set(2, 0, -0.775035); S_a.set(2, 1, -0.293982); S_a.set(2, 2, 0.559370);

  Matrix S_b;
  S_b.set(0, 0, -0.419979); S_b.set(0, 1, 0.301765); S_b.set(0, 2, -0.8558940);
  S_b.set(1, 0, -0.653011); S_b.set(1, 1, -0.755415); S_b.set(1, 2, 0.054087);
  S_b.set(2, 0, -0.630234); S_b.set(2, 1, 0.581624); S_b.set(2, 2, 0.514314);

  Matrix S_c = CGAL::Optimal_bounding_box::internal::mean(S_a, S_b, traits);
  check_equality(S_c(0,0), -0.35111);
  check_equality(S_c(0,1), 0.79308);
  check_equality(S_c(0,2), -0.49774);
  check_equality(S_c(1,0), -0.61398);
  check_equality(S_c(1,1), -0.59635);
  check_equality(S_c(1,2), -0.51710);
  check_equality(S_c(2,0), -0.70693);
  check_equality(S_c(2,1), 0.12405);
  check_equality(S_c(2,2), 0.69632);
}

void test_centroid(const Traits& traits)
{
  Matrix S_a;
  S_a.set(0, 0, -0.588443); S_a.set(0, 1, 0.807140); S_a.set(0, 2, -0.047542);
  S_a.set(1, 0, -0.786228); S_a.set(1, 1, -0.584933); S_a.set(1, 2, -0.199246);
  S_a.set(2, 0, -0.188629); S_a.set(2, 1, -0.079867); S_a.set(2, 2, 0.978795);

  Matrix S_b;
  S_b.set(0, 0, -0.2192721); S_b.set(0, 1, 0.2792986); S_b.set(0, 2, -0.9348326);
  S_b.set(1, 0, -0.7772152); S_b.set(1, 1, -0.6292092); S_b.set(1, 2, -0.005686);
  S_b.set(2, 0, -0.5897934); S_b.set(2, 1, 0.7253193); S_b.set(2, 2, 0.3550431);

  Matrix S_c;
  S_c.set(0, 0, -0.32657); S_c.set(0, 1, -0.60013); S_c.set(0, 2, -0.730206);
  S_c.set(1, 0, -0.20022); S_c.set(1, 1, -0.71110); S_c.set(1, 2, 0.67398);
  S_c.set(2, 0, -0.92372); S_c.set(2, 1, 0.36630); S_c.set(2, 2, 0.11207);

  Matrix S_centroid = CGAL::Optimal_bounding_box::internal::nm_centroid(S_a, S_b, S_c, traits);
  check_equality(S_centroid(0,0), -0.419979);
  check_equality(S_centroid(0,1), 0.301765);
  check_equality(S_centroid(0,2), -0.855894);
  check_equality(S_centroid(1,0), -0.653011);
  check_equality(S_centroid(1,1), -0.755415);
  check_equality(S_centroid(1,2), 0.054087);
  check_equality(S_centroid(2,0), -0.630234);
  check_equality(S_centroid(2,1), 0.581624);
  check_equality(S_centroid(2,2), 0.514314);
}

void test_nelder_mead(const Traits& traits)
{
  std::array<Point_3, 4> points;
  points[0] = Point_3(0.866802, 0.740808, 0.895304);
  points[1] = Point_3(0.912651, 0.761565, 0.160330);
  points[2] = Point_3(0.093661, 0.892578, 0.737412);
  points[3] = Point_3(0.166461, 0.149912, 0.364944);

  // one simplex
  std::array<Vertex, 4> simplex;

  Matrix v0, v1, v2, v3;
  v0.set(0, 0, -0.2192721); v0.set(0, 1, 0.2792986); v0.set(0, 2, -0.9348326);
  v0.set(1, 0, -0.7772152); v0.set(1, 1, -0.6292092); v0.set(1, 2, -0.0056861);
  v0.set(2, 0, -0.5897934); v0.set(2, 1, 0.7253193); v0.set(2, 2, 0.3550431);

  v1.set(0, 0, -0.588443); v1.set(0, 1, 0.807140); v1.set(0, 2, -0.047542);
  v1.set(1, 0, -0.786228); v1.set(1, 1, -0.584933); v1.set(1, 2, -0.199246);
  v1.set(2, 0, -0.188629); v1.set(2, 1, -0.079867); v1.set(2, 2, 0.978795);

  v2.set(0, 0, -0.277970); v2.set(0, 1, 0.953559); v2.set(0, 2, 0.116010);
  v2.set(1, 0, -0.567497); v2.set(1, 1, -0.065576); v2.set(1, 2, -0.820760);
  v2.set(2, 0, -0.775035); v2.set(2, 1, -0.293982); v2.set(2, 2, 0.559370);

  v3.set(0, 0, -0.32657); v3.set(0, 1, -0.60013); v3.set(0, 2, -0.73020);
  v3.set(1, 0, -0.20022); v3.set(1, 1, -0.71110); v3.set(1, 2, 0.67398);
  v3.set(2, 0, -0.92372); v3.set(2, 1, 0.36630); v3.set(2, 2, 0.11207);

  simplex[0] = Vertex{v0, points, traits};
  simplex[1] = Vertex{v1, points, traits};
  simplex[2] = Vertex{v2, points, traits};
  simplex[3] = Vertex{v3, points, traits};

  std::size_t nm_iterations = 19;
  CGAL::Optimal_bounding_box::internal::nelder_mead(simplex, nm_iterations, points, traits);

  const Matrix& v0_new = simplex[0].matrix();
  check_equality(v0_new(0,0), -0.288975);
  check_equality(v0_new(0,1), 0.7897657);
  check_equality(v0_new(0,2), -0.541076);
  check_equality(v0_new(1,0), -0.9407046);
  check_equality(v0_new(1,1), -0.3391466);
  check_equality(v0_new(1,2), 0.0073817);
  check_equality(v0_new(2,0), -0.1776743);
  check_equality(v0_new(2,1), 0.5111260);
  check_equality(v0_new(2,2), 0.84094);

  const Matrix& v1_new = simplex[1].matrix();
  check_equality(v1_new(0,0), -0.458749);
  check_equality(v1_new(0,1), 0.823283);
  check_equality(v1_new(0,2), -0.334296);
  check_equality(v1_new(1,0), -0.885235);
  check_equality(v1_new(1,1), -0.455997);
  check_equality(v1_new(1,2), 0.091794);
  check_equality(v1_new(2,0), -0.076866);
  check_equality(v1_new(2,1), 0.338040);
  check_equality(v1_new(2,2), 0.937987);

  const Matrix& v2_new = simplex[2].matrix();
  check_equality(v2_new(0,0), -0.346582);
  check_equality(v2_new(0,1), 0.878534);
  check_equality(v2_new(0,2), -0.328724);
  check_equality(v2_new(1,0), -0.936885);
  check_equality(v2_new(1,1), -0.341445);
  check_equality(v2_new(1,2), 0.075251);
  check_equality(v2_new(2,0), -0.046131);
  check_equality(v2_new(2,1), 0.334057);
  check_equality(v2_new(2,2), 0.941423);

  const Matrix& v3_new = simplex[3].matrix();
  check_equality(v3_new(0,0), -0.394713);
  check_equality(v3_new(0,1), 0.791782);
  check_equality(v3_new(0,2), -0.466136);
  check_equality(v3_new(1,0), -0.912112);
  check_equality(v3_new(1,1), -0.398788);
  check_equality(v3_new(1,2), 0.094972);
  check_equality(v3_new(2,0), -0.110692);
  check_equality(v3_new(2,1), 0.462655);
  check_equality(v3_new(2,2), 0.879601);
}

int main(int, char**)
{
  Traits traits;

  test_simplex_operations(traits);
  test_centroid(traits);
  test_nelder_mead(traits);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
