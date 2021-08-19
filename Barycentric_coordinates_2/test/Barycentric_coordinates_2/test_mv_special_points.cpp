#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>
#include <boost/math/special_functions/fpclassify.hpp>

template<typename Kernel>
void test_mv_special_points() {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 1),
    Point_2(FT(7) / FT(4), FT(1) / FT(2)),
    Point_2(FT(7) / FT(4), FT(5) / FT(2)),
    Point_2(1, 2),
    Point_2(0, 3),
    Point_2(FT(1) / FT(2), FT(3) / FT(2))
  };

  const Point_2 queries[11] = {
    Point_2(FT(1) + (FT(1) / FT(std::pow(10.0, 300.0))), FT(2) - (FT(1) / FT(std::pow(10.0, 300.0)))),
    Point_2(FT(1) + (FT(1) / FT(std::pow(10.0, 300.0))), FT(1) + (FT(1) / FT(std::pow(10.0, 300.0)))),
    Point_2(1, FT(3) / FT(2)),
    Point_2(FT(5) / FT(4), FT(5) / FT(4)),
    Point_2(FT(5) / FT(4), FT(7) / FT(4)),
    Point_2(FT(3) / FT(2), FT(3) / FT(2)),

    Point_2(FT(7) / FT(4) - (FT(1) / FT(std::pow(10.0, 300.0))), FT(7) / FT(4) - (FT(1) / FT(std::pow(10.0, 300.0)))),
    Point_2(FT(7) / FT(4) - (FT(1) / FT(std::pow(10.0, 300.0))), FT(5) / FT(4) + (FT(1) / FT(std::pow(10.0, 300.0)))),

    Point_2(FT(3) / FT(4) - (FT(1) / FT(std::pow(10.0, 300.0))), FT(3) / FT(4) + (FT(1) / FT(std::pow(10.0, 300.0)))),
    Point_2(FT(3) / FT(4) - (FT(1) / FT(std::pow(10.0, 300.0))), FT(9) / FT(4) - (FT(1) / FT(std::pow(10.0, 300.0)))),
    Point_2(FT(1) / FT(2) + (FT(1) / FT(std::pow(10.0, 300.0))), FT(3) / FT(2))
  };

  std::size_t count = 0;
  const FT epsilon = FT(1) / FT(1000000000000000);

  std::vector<FT> coordinates;
  for (std::size_t i = 0; i < 11; ++i) {
    CGAL::Barycentric_coordinates::mean_value_coordinates_2(
      vertices, queries[i], std::back_inserter(coordinates));

    assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 0])));
    assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 0])));

    assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 1])));
    assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 1])));

    assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 2])));
    assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 2])));

    assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 3])));
    assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 3])));

    assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 4])));
    assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 4])));

    assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 5])));
    assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 5])));

    assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 6])));
    assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 6])));

    const FT coordinate_sum =
      coordinates[count + 0] +
      coordinates[count + 1] +
      coordinates[count + 2] +
      coordinates[count + 3] +
      coordinates[count + 4] +
      coordinates[count + 5] +
      coordinates[count + 6] ;
    const Point_2 linear_combination(
      vertices[0].x() * coordinates[count + 0] +
      vertices[1].x() * coordinates[count + 1] +
      vertices[2].x() * coordinates[count + 2] +
      vertices[3].x() * coordinates[count + 3] +
      vertices[4].x() * coordinates[count + 4] +
      vertices[5].x() * coordinates[count + 5] +
      vertices[6].x() * coordinates[count + 6] ,
      vertices[0].y() * coordinates[count + 0] +
      vertices[1].y() * coordinates[count + 1] +
      vertices[2].y() * coordinates[count + 2] +
      vertices[3].y() * coordinates[count + 3] +
      vertices[4].y() * coordinates[count + 4] +
      vertices[5].y() * coordinates[count + 5] +
      vertices[6].y() * coordinates[count + 6] );
    const Point_2 difference(
      linear_combination.x() - queries[i].x(),
      linear_combination.y() - queries[i].y());
    assert(
      CGAL::abs(coordinate_sum - FT(1)) < epsilon &&
      CGAL::abs(difference.x()) < epsilon &&
      CGAL::abs(difference.y()) < epsilon );
    count += 7;
  }
}

int main() {

  test_mv_special_points< CGAL::Simple_cartesian<double> >();
  test_mv_special_points< CGAL::Exact_predicates_inexact_constructions_kernel >();

  std::cout << "test_mv_special_points: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
