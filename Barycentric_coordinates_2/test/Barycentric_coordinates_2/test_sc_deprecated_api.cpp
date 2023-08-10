#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/segment_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Segment_coordinates_2<Kernel> Segment_coordinates;

typedef std::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
  const Point first_vertex  = Point(0, 0);
  const Point second_vertex = Point(1, 0);

  Segment_coordinates segment_coordinates(first_vertex, second_vertex);

  const Point query_points[6] = {
    Point(Scalar(2) /Scalar(5) , 0),
    Point(1                    , 0),
    Point(Scalar(7) /Scalar(10), 0),
    Point(Scalar(-3)/Scalar(10), 0),
    Point(Scalar(6) /Scalar(5) , 0),
    Point(0                    , 0) };

  const Scalar expected_coordinates[12] = {
    Scalar(3) /Scalar(5) , Scalar(2) /Scalar(5) ,
    0                    , 1                    ,
    Scalar(3) /Scalar(10), Scalar(7) /Scalar(10),
    Scalar(13)/Scalar(10), Scalar(-3)/Scalar(10),
    Scalar(-1)/Scalar(5) , Scalar(6) /Scalar(5) ,
    1                    , 0                     };

  Coordinate_vector old_coordinates;
  Coordinate_vector new_coordinates;

  int count = 0;
  for (int i = 0; i < 6; ++i) {
    segment_coordinates(query_points[i], std::back_inserter(old_coordinates));
    CGAL::Barycentric_coordinates::segment_coordinates_2(
      first_vertex, second_vertex, query_points[i], std::back_inserter(new_coordinates));

    assert(
      old_coordinates[count + 0] - expected_coordinates[count + 0] == Scalar(0) &&
      old_coordinates[count + 1] - expected_coordinates[count + 1] == Scalar(0) );

    assert(
      old_coordinates[count + 0] - new_coordinates[count + 0] == Scalar(0) &&
      old_coordinates[count + 1] - new_coordinates[count + 1] == Scalar(0) );

    if (
      old_coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
      old_coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0)  )
    {
      cout << endl << "Segment_coordinates_deprecated_api_test: FAILED." << endl << endl;
      exit(EXIT_FAILURE);
    }
    count += 2;
  }
  old_coordinates.clear();

  count = 0;
  for (int i = 0; i < 6; ++i) {
    const auto bc = CGAL::Barycentric_coordinates::
      compute_segment_coordinates_2(first_vertex, second_vertex, query_points[i], Kernel());

    assert(
      bc[0] - expected_coordinates[count + 0] == Scalar(0) &&
      bc[1] - expected_coordinates[count + 1] == Scalar(0) );

    if (
      bc[0] - expected_coordinates[count + 0] != Scalar(0) ||
      bc[1] - expected_coordinates[count + 1] != Scalar(0)  )
    {
      cout << endl << "Segment_coordinates_deprecated_api_test: FAILED." << endl << endl;
      exit(EXIT_FAILURE);
    }
    count += 2;
  }

  cout << endl << "Segment_coordinates_deprecated_api_test: PASSED." << endl << endl;
  return EXIT_SUCCESS;
}
