// Author(s) : Dmitry Anisimov.
// In this test we compute Discrete Harmonic weights for ~9800 strictly interior points with respect to
// a pentagon. Then we sum them up and normalize by this sum. What we expect is Discrete Harmonic coordinates.
// The chosen data type is exact. Some points very close to the boundary are also used.

// Does not work with inexact kernel. Because weights cannot be computed with a distance 1.0e-300 away from the boundary.

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
    const Point vertices[5] = { Point(0                   , 0                  ),
                                Point(1                   , 0                  ),
                                Point(Scalar(5) /Scalar(4), Scalar(3)/Scalar(4)),
                                Point(Scalar(1) /Scalar(2), Scalar(3)/Scalar(2)),
                                Point(Scalar(-1)/Scalar(4), Scalar(3)/Scalar(4))
                              };

    const Polygon pentagon(vertices, vertices + 5);

    Discrete_harmonic_coordinates discrete_harmonic_coordinates(pentagon);

    Coordinate_vector weights;
    Coordinate_vector coordinates;
    Coordinate_vector expected_coordinates;

    const Scalar step  = Scalar(1)  / Scalar(100);
    const Scalar scale = Scalar(100);

    int count = 0;
    const Scalar limit = scale*step;

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(x, y);

            const Output_type w_result = discrete_harmonic_coordinates.compute_weights(point, weights);

            Scalar W = Scalar(0);
            for(int j = 0; j < 5; ++j) W += weights[count + j];

            const Scalar inverted_W = Scalar(1) / W;

            for(int j = 0; j < 5; ++j) coordinates.push_back(weights[count + j] * inverted_W);

            const Output_type c_result = discrete_harmonic_coordinates.compute(point, expected_coordinates);

            if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
                coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
                coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0) ||
                coordinates[count + 3] - expected_coordinates[count + 3] != Scalar(0) ||
                coordinates[count + 4] - expected_coordinates[count + 4] != Scalar(0)  )
            {
                cout << endl << w_result.second << " " << c_result.second << " DH_weights_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 5;
        }
    }

    const Point query_points[5] = { Point(Scalar(1) / Scalar(pow(10.0, 300.0))          , Scalar(1)/Scalar(pow(10.0, 300.0))),
                                    Point(Scalar(1) / Scalar(4)                         , Scalar(1)/Scalar(pow(10.0, 300.0))),
                                    Point(Scalar(1) / Scalar(2)                         , Scalar(1)/Scalar(pow(10.0, 300.0))),
                                    Point(Scalar(3) / Scalar(4)                         , Scalar(1)/Scalar(pow(10.0, 300.0))),
                                    Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 300.0)))
                                  };

    for(int i = 0; i < 5; ++i) {
        const Output_type w_result = discrete_harmonic_coordinates.compute_weights(query_points[i], weights);

        Scalar W = Scalar(0);
        for(int j = 0; j < 5; ++j) W += weights[count + j];

        const Scalar inverted_W = Scalar(1) / W;

        for(int j = 0; j < 5; ++j) coordinates.push_back(weights[count + j] * inverted_W);

        const Output_type c_result = discrete_harmonic_coordinates.compute(query_points[i], expected_coordinates);

        if( coordinates[count + 0] - expected_coordinates[count + 0] != Scalar(0) ||
            coordinates[count + 1] - expected_coordinates[count + 1] != Scalar(0) ||
            coordinates[count + 2] - expected_coordinates[count + 2] != Scalar(0) ||
            coordinates[count + 3] - expected_coordinates[count + 3] != Scalar(0) ||
            coordinates[count + 4] - expected_coordinates[count + 4] != Scalar(0)  )
        {
            cout << endl << w_result.second << " " << c_result.second << " DH_weights_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 5;
    }

    cout << endl << "DH_weights_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
