// Author(s) : Dmitry Anisimov.
// We use a simple concave polygon and an exact data type
// in order to test coordinates computed for points along the boundary of the polygon.
// As a base coordinate function, we take Mean Value coordinates.
// We test all compute_at_vertex(), compute_on_edge(), and compute() functions with std::back_inserter(std::vector) output.

// Works with inexact kernel, too.

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Mean_value_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::MV_coordinates_2<Polygon, Coordinate_vector> Mean_value_coordinates;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point vertices[7] = { Point(0                   , 0                  ),
                                Point(1                   , 0                  ),
                                Point(Scalar(1) /Scalar(2), 1                  ),
                                Point(Scalar(3) /Scalar(2), Scalar(3)/Scalar(2)),
                                Point(Scalar(-1)/Scalar(2), Scalar(3)/Scalar(2)),
                                Point(0                   , 1                  ),
                                Point(Scalar(-1)/Scalar(2), Scalar(1)/Scalar(2))
                              };

    const Polygon concave_polygon(vertices, vertices + 7);

    Coordinate_vector coordinates;

    Mean_value_coordinates mean_value_coordinates(concave_polygon);

    const Point query_points[7] = { Point(Scalar(1) /Scalar(2), 0                  ),
                                    Point(Scalar(3) /Scalar(4), Scalar(1)/Scalar(2)),
                                    Point(Scalar(1) /Scalar(2), 1                  ),
                                    Point(Scalar(3) /Scalar(4), Scalar(3)/Scalar(2)),
                                    Point(0                   , 1                  ),
                                    Point(Scalar(-1)/Scalar(8), Scalar(7)/Scalar(8)),
                                    Point(Scalar(-1)/Scalar(8), Scalar(1)/Scalar(8))
                                  };

    const Scalar expected_coordinates[49] = { Scalar(1)/Scalar(2), Scalar(1)/Scalar(2), 0, 0, 0, 0, 0,
                                              0, Scalar(1)/Scalar(2), Scalar(1)/Scalar(2), 0, 0, 0, 0,
                                              0, 0                  , 1                  , 0, 0, 0, 0,
                                              0, 0, 0, Scalar(5)/Scalar(8), Scalar(3)/Scalar(8), 0, 0,
                                              0, 0, 0, 0, 0                                    , 1, 0,
                                              0, 0, 0, 0, 0, Scalar(3)/Scalar(4), Scalar(1)/Scalar(4),
                                              Scalar(3)/Scalar(4), 0, 0, 0, 0, 0, Scalar(1)/Scalar(4)
                                            };

    int count = 0;
    Output_type result = mean_value_coordinates.compute_at_vertex(2, coordinates);

    if( coordinates[count + 0] - expected_coordinates[14 + 0] != Scalar(0) ||
        coordinates[count + 1] - expected_coordinates[14 + 1] != Scalar(0) ||
        coordinates[count + 2] - expected_coordinates[14 + 2] != Scalar(0) ||
        coordinates[count + 3] - expected_coordinates[14 + 3] != Scalar(0) ||
        coordinates[count + 4] - expected_coordinates[14 + 4] != Scalar(0) ||
        coordinates[count + 5] - expected_coordinates[14 + 5] != Scalar(0) ||
        coordinates[count + 6] - expected_coordinates[14 + 6] != Scalar(0)  )
    {
        cout << endl << result.second << " Computation_on_boundary_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }
    count += 7;

    result = mean_value_coordinates.compute_at_vertex(5, coordinates);

    if( coordinates[count + 0] - expected_coordinates[28 + 0] != Scalar(0) ||
        coordinates[count + 1] - expected_coordinates[28 + 1] != Scalar(0) ||
        coordinates[count + 2] - expected_coordinates[28 + 2] != Scalar(0) ||
        coordinates[count + 3] - expected_coordinates[28 + 3] != Scalar(0) ||
        coordinates[count + 4] - expected_coordinates[28 + 4] != Scalar(0) ||
        coordinates[count + 5] - expected_coordinates[28 + 5] != Scalar(0) ||
        coordinates[count + 6] - expected_coordinates[28 + 6] != Scalar(0)  )
    {
        cout << endl << result.second << " Computation_on_boundary_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }
    count += 7;

    result = mean_value_coordinates.compute(query_points[2], coordinates, CGAL::Barycentric_coordinates::AT_VERTEX);

    if( coordinates[count + 0] - expected_coordinates[14 + 0] != Scalar(0) ||
        coordinates[count + 1] - expected_coordinates[14 + 1] != Scalar(0) ||
        coordinates[count + 2] - expected_coordinates[14 + 2] != Scalar(0) ||
        coordinates[count + 3] - expected_coordinates[14 + 3] != Scalar(0) ||
        coordinates[count + 4] - expected_coordinates[14 + 4] != Scalar(0) ||
        coordinates[count + 5] - expected_coordinates[14 + 5] != Scalar(0) ||
        coordinates[count + 6] - expected_coordinates[14 + 6] != Scalar(0)  )
    {
        cout << endl << result.second << " Computation_on_boundary_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }
    count += 7;

    result = mean_value_coordinates.compute(query_points[4], coordinates, CGAL::Barycentric_coordinates::AT_VERTEX);

    if( coordinates[count + 0] - expected_coordinates[28 + 0] != Scalar(0) ||
        coordinates[count + 1] - expected_coordinates[28 + 1] != Scalar(0) ||
        coordinates[count + 2] - expected_coordinates[28 + 2] != Scalar(0) ||
        coordinates[count + 3] - expected_coordinates[28 + 3] != Scalar(0) ||
        coordinates[count + 4] - expected_coordinates[28 + 4] != Scalar(0) ||
        coordinates[count + 5] - expected_coordinates[28 + 5] != Scalar(0) ||
        coordinates[count + 6] - expected_coordinates[28 + 6] != Scalar(0)  )
    {
        cout << endl << result.second << " Computation_on_boundary_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    coordinates.clear();

    count = 0;
    for(int i = 0; i < 7; ++i) {
        result = mean_value_coordinates.compute_on_edge(query_points[i], i, coordinates);

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0) ||
            coordinates[count + 3] - expected_coordinates[count + 3] != Scalar(0) ||
            coordinates[count + 4] - expected_coordinates[count + 4] != Scalar(0) ||
            coordinates[count + 5] - expected_coordinates[count + 5] != Scalar(0) ||
            coordinates[count + 6] - expected_coordinates[count + 6] != Scalar(0)  )
        {
            cout << endl << result.second << " Computation_on_boundary_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 7;
    }
    coordinates.clear();

    count = 0;
    for(int i = 0; i < 7; ++i) {
        result = mean_value_coordinates.compute(query_points[i], coordinates, CGAL::Barycentric_coordinates::ON_BOUNDARY);

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0) ||
            coordinates[count + 3] - expected_coordinates[count + 3] != Scalar(0) ||
            coordinates[count + 4] - expected_coordinates[count + 4] != Scalar(0) ||
            coordinates[count + 5] - expected_coordinates[count + 5] != Scalar(0) ||
            coordinates[count + 6] - expected_coordinates[count + 6] != Scalar(0)  )
        {
            cout << endl << result.second << " Computation_on_boundary_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 7;
    }

    cout << endl << "Computation_on_boundary_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
