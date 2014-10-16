// Author(s) : Dmitry Anisimov.
// We use a simple non-regular strictly convex hexagon and an exact data type
// in order to test coordinates computed for points along all the open edges of the polygon.
// As a base coordinate function, we take Discrete Harmonic coordinates.
// We test both compute_on_edge() and compute() functions with std::vector output.

// Works with inexact kernel, too.

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Discrete_harmonic_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::DH_coordinates_2<Polygon, Coordinate_vector> Discrete_harmonic_coordinates;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point vertices[6] = { Point(0                   , 0                  ),
                                Point(1                   , 0                  ),
                                Point(Scalar(3) /Scalar(2), 1                  ),
                                Point(Scalar(1) /Scalar(2), 2                  ),
                                Point(Scalar(-1)/Scalar(2), Scalar(3)/Scalar(2)),
                                Point(Scalar(-1)/Scalar(2), Scalar(1)/Scalar(2))
                              };

    const Polygon hexagon(vertices, vertices + 6);

    Coordinate_vector coordinates;

    Discrete_harmonic_coordinates discrete_harmonic_coordinates(hexagon);

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

    int count = 0;
    for(int i = 0; i < 6; ++i) {
        const Output_type result = discrete_harmonic_coordinates.compute_on_edge(query_points[i], i, coordinates);

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0) ||
            coordinates[count + 3] - expected_coordinates[count + 3] != Scalar(0) ||
            coordinates[count + 4] - expected_coordinates[count + 4] != Scalar(0) ||
            coordinates[count + 5] - expected_coordinates[count + 5] != Scalar(0)  )
        {
            cout << endl << result.second << " Computation_on_edges_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 6;
    }
    coordinates.clear();

    count = 0;
    for(int i = 0; i < 6; ++i) {
        const Output_type result = discrete_harmonic_coordinates.compute(query_points[i], coordinates, CGAL::Barycentric_coordinates::ON_BOUNDARY);

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0) ||
            coordinates[count + 3] - expected_coordinates[count + 3] != Scalar(0) ||
            coordinates[count + 4] - expected_coordinates[count + 4] != Scalar(0) ||
            coordinates[count + 5] - expected_coordinates[count + 5] != Scalar(0)  )
        {
            cout << endl << result.second << " Computation_on_edges_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 6;
    }

    cout << endl << "Computation_on_edges_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
