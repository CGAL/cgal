#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::vector<Point>  Point_vector;

typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;
typedef CGAL::Barycentric_coordinates::Mean_value_2<Kernel> Mean_value;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel> Mean_value_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
  const Point first_vertex  = Point(0, 0);
  const Point second_vertex = Point(1, 0);
  const Point third_vertex  = Point(0, 1);

  Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);

  Point_vector vertices(3);

  vertices[0] = Point(0, 0); vertices[1] = Point(1, 0); vertices[2] = Point(0, 1);

  Mean_value_coordinates mean_value_coordinates(vertices.begin(), vertices.end());

  Coordinate_vector tri_coordinates;
  Coordinate_vector old_coordinates;
  Coordinate_vector new_coordinates;

  const Scalar step  = Scalar(1) / Scalar(100);
  const Scalar scale = Scalar(50);
  const Scalar tol   = Scalar(1) / Scalar(10000000000);

  int count = 0;
  const Scalar limit = scale * step;
  const Scalar epsilon = Scalar(1) / Scalar(std::pow(10.0, 14.0));

  for (Scalar x = step; x < limit; x += step) {
    for (Scalar y = step; y < limit; y += step) {
      const Point point(x, y);

      const Output_type tri_result = triangle_coordinates(point, tri_coordinates);
      const Output_type  mv_result = mean_value_coordinates(point, old_coordinates);
      CGAL::Barycentric_coordinates::mean_value_coordinates_2(
        vertices, point, std::back_inserter(new_coordinates));

      assert(
        (tri_coordinates[count + 0] - old_coordinates[count + 0]) < epsilon &&
        (tri_coordinates[count + 1] - old_coordinates[count + 1]) < epsilon &&
        (tri_coordinates[count + 2] - old_coordinates[count + 2]) < epsilon );

      assert(
        CGAL::abs(old_coordinates[count + 0] - new_coordinates[count + 0]) < tol &&
        CGAL::abs(old_coordinates[count + 1] - new_coordinates[count + 1]) < tol &&
        CGAL::abs(old_coordinates[count + 2] - new_coordinates[count + 2]) < tol );

      if (
        (tri_coordinates[count + 0] - old_coordinates[count + 0]) > epsilon ||
        (tri_coordinates[count + 1] - old_coordinates[count + 1]) > epsilon ||
        (tri_coordinates[count + 2] - old_coordinates[count + 2]) > epsilon  )
      {
        cout << endl << "MV_deprecated_api_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
      }
      count += 3;
    }
  }

  cout << endl << "MV_deprecated_api_test: PASSED." << endl << endl;
  return EXIT_SUCCESS;
}
