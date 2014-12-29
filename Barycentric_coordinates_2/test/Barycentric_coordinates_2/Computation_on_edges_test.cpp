// Author: Dmitry Anisimov.
// We use a simple non-regular strictly convex hexagon and an exact data type
// in order to test coordinates computed for points along all the open edges of the polygon.
// As a base coordinate function, we take discrete harmonic coordinates.
// We test both compute_on_edge() and operator() functions with std::vector output.

// Works with inexact kernel, too.

#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::vector<Point>  Point_vector;

typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Discrete_harmonic_2<Kernel> Discrete_harmonic;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Discrete_harmonic, Kernel> Discrete_harmonic_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    Point_vector vertices(6);

    vertices[0] = Point(0, 0);                   vertices[1] = Point(1, 0);                                      vertices[2] = Point(Scalar(3) /Scalar(2), 1);
    vertices[3] = Point(Scalar(1)/Scalar(2), 2); vertices[4] = Point(Scalar(-1)/Scalar(2), Scalar(3)/Scalar(2)); vertices[5] = Point(Scalar(-1)/Scalar(2), Scalar(1)/Scalar(2));

    Discrete_harmonic_coordinates discrete_harmonic_coordinates(vertices.begin(), vertices.end());

    const Point query_points[6] = { Point(Scalar(1) /Scalar(2), 0                   ),
                                    Point(Scalar(5) /Scalar(4), Scalar(1) /Scalar(2)),
                                    Point(Scalar(3) /Scalar(4), Scalar(7) /Scalar(4)),
                                    Point(0                   , Scalar(7) /Scalar(4)),
                                    Point(Scalar(-1)/Scalar(2), Scalar(11)/Scalar(8)),
                                    Point(Scalar(-1)/Scalar(4), Scalar(1) /Scalar(4))
                                  };

    const Scalar expected_coordinates[36] = { Scalar(1)/Scalar(2), Scalar(1)/Scalar(2), 0, 0, 0, 0,
                                              0, Scalar(1)/Scalar(2), Scalar(1)/Scalar(2), 0, 0, 0,
                                              0, 0, Scalar(1)/Scalar(4), Scalar(3)/Scalar(4), 0, 0,
                                              0, 0, 0, Scalar(1)/Scalar(2), Scalar(1)/Scalar(2), 0,
                                              0, 0, 0, 0, Scalar(7)/Scalar(8), Scalar(1)/Scalar(8),
                                              Scalar(1)/Scalar(2), 0, 0, 0, 0, Scalar(1)/Scalar(2)
                                            };

    Coordinate_vector coordinates;

    int count = 0;
    for(int i = 0; i < 6; ++i) {
        const Output_type result = discrete_harmonic_coordinates.compute_on_edge(query_points[i], i, coordinates);

        assert(coordinates[count + 0] - expected_coordinates[count + 0] == Scalar(0) && coordinates[count + 1] - expected_coordinates[count + 1] == Scalar(0) &&
               coordinates[count + 2] - expected_coordinates[count + 2] == Scalar(0) && coordinates[count + 3] - expected_coordinates[count + 3] == Scalar(0) &&
               coordinates[count + 4] - expected_coordinates[count + 4] == Scalar(0) && coordinates[count + 5] - expected_coordinates[count + 5] == Scalar(0) );

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0) ||
            coordinates[count + 3] - expected_coordinates[count + 3] != Scalar(0) ||
            coordinates[count + 4] - expected_coordinates[count + 4] != Scalar(0) ||
            coordinates[count + 5] - expected_coordinates[count + 5] != Scalar(0)  )
        {
            cout << endl << "Computation_on_edges_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 6;
    }
    coordinates.clear();

    count = 0;
    for(int i = 0; i < 6; ++i) {
        const Output_type result = discrete_harmonic_coordinates(query_points[i], coordinates, CGAL::Barycentric_coordinates::ON_BOUNDARY);

        assert(coordinates[count + 0] - expected_coordinates[count + 0] == Scalar(0) && coordinates[count + 1] - expected_coordinates[count + 1] == Scalar(0) &&
               coordinates[count + 2] - expected_coordinates[count + 2] == Scalar(0) && coordinates[count + 3] - expected_coordinates[count + 3] == Scalar(0) &&
               coordinates[count + 4] - expected_coordinates[count + 4] == Scalar(0) && coordinates[count + 5] - expected_coordinates[count + 5] == Scalar(0) );

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0) ||
            coordinates[count + 3] - expected_coordinates[count + 3] != Scalar(0) ||
            coordinates[count + 4] - expected_coordinates[count + 4] != Scalar(0) ||
            coordinates[count + 5] - expected_coordinates[count + 5] != Scalar(0)  )
        {
            cout << endl << "Computation_on_edges_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 6;
    }

    cout << endl << "Computation_on_edges_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
