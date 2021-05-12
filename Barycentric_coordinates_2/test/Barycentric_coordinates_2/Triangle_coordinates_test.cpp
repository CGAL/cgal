// Author: Dmitry Anisimov.
// We use a simple triangle and an exact data type
// in order to test coordinates computed for some points in the plane (inside triangle and outside it).
// We test the function operator() and its overload.

// Works with inexact kernel, too.

#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point first_vertex  = Point(0, 0);
    const Point second_vertex = Point(1, 0);
    const Point third_vertex  = Point(Scalar(1)/Scalar(2), 1);

    Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);

    const Point query_points[5] = { Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(2)),
                                    Point(Scalar(1)/Scalar(2), 0                  ),
                                    Point(Scalar(1)/Scalar(2), 1                  ),
                                    Point(1                  , Scalar(1)/Scalar(2)),
                                    Point(0                  , Scalar(1)/Scalar(2))
                                  };

    const Scalar expected_coordinates[15] = { Scalar(1) /Scalar(4), Scalar(1) /Scalar(4), Scalar(1)/Scalar(2),
                                              Scalar(1) /Scalar(2), Scalar(1) /Scalar(2), 0                  ,
                                              0                   , 0                   , 1                  ,
                                              Scalar(-1)/Scalar(4), Scalar(3) /Scalar(4), Scalar(1)/Scalar(2),
                                              Scalar(3) /Scalar(4), Scalar(-1)/Scalar(4), Scalar(1)/Scalar(2)
                                            };

    Coordinate_vector coordinates;

    int count = 0;
    for(int i = 0; i < 5; ++i) {
        const Output_type result = triangle_coordinates(query_points[i], std::back_inserter(coordinates));

        assert(coordinates[count + 0] - expected_coordinates[count + 0] == Scalar(0) &&
               coordinates[count + 1] - expected_coordinates[count + 1] == Scalar(0) &&
               coordinates[count + 2] - expected_coordinates[count + 2] == Scalar(0) );

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0)  )
        {
            cout << endl << "Triangle_coordinates_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 3;
    }
    coordinates.clear();

    count = 0;
    for(int i = 0; i < 5; ++i) {
        const Output_type result = triangle_coordinates(query_points[i], coordinates);

        assert(coordinates[count + 0] - expected_coordinates[count + 0] == Scalar(0) &&
               coordinates[count + 1] - expected_coordinates[count + 1] == Scalar(0) &&
               coordinates[count + 2] - expected_coordinates[count + 2] == Scalar(0) );

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0)  )
        {
            cout << endl << "Triangle_coordinates_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 3;
    }

    cout << endl << "Triangle_coordinates_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
