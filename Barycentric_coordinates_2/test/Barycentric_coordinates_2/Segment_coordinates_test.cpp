// Author: Dmitry Anisimov.
// We use a simple segment of length 1 and an exact data type
// in order to test coordinates computed for some points along the segment (inside and outside).
// We test the function operator() and its overload.

// Does not work with inexact kernel. Get inconsistency when comparing coordinates with expected_coordinates.

#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT        Scalar;
typedef Kernel::Point_2   Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Segment_coordinates_2<Kernel> Segment_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point first_vertex  = Point(0, 0);
    const Point second_vertex = Point(1, 0);

    Segment_coordinates segment_coordinates(first_vertex, second_vertex);

    const Point query_points[6] = { Point(Scalar(2) /Scalar(5) , 0),
                                    Point(1                    , 0),
                                    Point(Scalar(7) /Scalar(10), 0),
                                    Point(Scalar(-3)/Scalar(10), 0),
                                    Point(Scalar(6) /Scalar(5) , 0),
                                    Point(0                    , 0)
                                  };

    const Scalar expected_coordinates[12] = { Scalar(3) /Scalar(5) , Scalar(2) /Scalar(5) ,
                                              0                    , 1                    ,
                                              Scalar(3) /Scalar(10), Scalar(7) /Scalar(10),
                                              Scalar(13)/Scalar(10), Scalar(-3)/Scalar(10),
                                              Scalar(-1)/Scalar(5) , Scalar(6) /Scalar(5) ,
                                              1                    , 0
                                            };

    Coordinate_vector coordinates;

    int count = 0;
    for(int i = 0; i < 6; ++i) {
        const Output_type result = segment_coordinates(query_points[i], std::back_inserter(coordinates));

        assert(coordinates[count + 0] - expected_coordinates[count + 0] == Scalar(0) &&
               coordinates[count + 1] - expected_coordinates[count + 1] == Scalar(0) );

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0)  )
        {
            cout << endl << "Segment_coordinates_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 2;
    }
    coordinates.clear();

    count = 0;
    for(int i = 0; i < 6; ++i) {
        const Output_type result = segment_coordinates(query_points[i], coordinates);

        assert(coordinates[count + 0] - expected_coordinates[count + 0] == Scalar(0) &&
               coordinates[count + 1] - expected_coordinates[count + 1] == Scalar(0) );

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0)  )
        {
            cout << endl << "Segment_coordinates_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 2;
    }

    cout << endl << "Segment_coordinates_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
