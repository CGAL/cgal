// Author(s) : Dmitry Anisimov.
// We use a simple triangle and an exact data type
// in order to test coordinates computed for some points in the plane (inside triangle and outside it).
// We test the function compute() and its overload.

// Works with inexact kernel, too.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Triangle_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;
typedef Kernel::Triangle_2 Triangle;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Triangle, Vector_insert_iterator> Triangle_coordinates;
typedef CGAL::Barycentric_coordinates::Tri_coordinates_2<Triangle, Coordinate_vector>  Triangle_coordinates_overload;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Triangle triangle( Point(0, 0), Point(1, 0), Point(Scalar(1)/Scalar(2), 1) );

    Coordinate_vector coordinates;

    Triangle_coordinates triangle_coordinates(triangle);
    Triangle_coordinates_overload triangle_coordinates_overload(triangle);

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

    int count = 0;
    for(int i = 0; i < 5; ++i) {
        const Output_type result = triangle_coordinates.compute(query_points[i], std::back_inserter(coordinates));

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0)  )
        {
            cout << endl << result.second << " Triangle_coordinates_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 3;
    }
    coordinates.clear();

    count = 0;
    for(int i = 0; i < 5; ++i) {
        const Output_type result = triangle_coordinates_overload.compute(query_points[i], coordinates);

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0)  )
        {
            cout << endl << result.second << " Triangle_coordinates_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 3;
    }

    cout << endl << "Triangle_coordinates_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
