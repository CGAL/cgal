// Author: Dmitry Anisimov.
// In this test we compute mean value weights for ~9800 strictly interior points with respect to
// a pentagon. Then we sum them up and normalize by this sum. What we expect is mean value coordinates.
// The chosen data type is exact. The used epsilon is 1.0e-14.

// Does not work with inexact kernel. We get inconsistency when comparing coordinates and expected_coordinates.

#include <cmath>
#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::vector<Point>  Point_vector;

typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Mean_value_2<Kernel> Mean_value;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel> Mean_value_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    Point_vector vertices(5);

    vertices[0] = Point(0, 0);                                      vertices[1] = Point(1, 0);                                      
    vertices[2] = Point(Scalar(5) /Scalar(4), Scalar(3)/Scalar(4)); vertices[3] = Point(Scalar(1)/Scalar(2), Scalar(3)/Scalar(2)); 
    vertices[4] = Point(Scalar(-1)/Scalar(4), Scalar(3)/Scalar(4));

    Mean_value_coordinates mean_value_coordinates(vertices.begin(), vertices.end());

    Coordinate_vector weights;
    Coordinate_vector coordinates;
    Coordinate_vector expected_coordinates;

    const Scalar step  = Scalar(1)  / Scalar(100);
    const Scalar scale = Scalar(100);

    const Scalar epsilon = Scalar(1) / Scalar(std::pow(10.0, 14.0));

    int count = 0;
    const Scalar limit = scale*step;

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(x, y);

            const Output_type w_result = mean_value_coordinates.compute_weights(point, weights);

            Scalar W = Scalar(0);
            for(int j = 0; j < 5; ++j) W += weights[count + j];

            const Scalar inverted_W = Scalar(1) / W;

            for(int j = 0; j < 5; ++j) coordinates.push_back(weights[count + j] * inverted_W);

            const Output_type c_result = mean_value_coordinates(point, expected_coordinates);

            assert(coordinates[count + 0] - expected_coordinates[count + 0] < epsilon &&
                   coordinates[count + 1] - expected_coordinates[count + 1] < epsilon &&
                   coordinates[count + 2] - expected_coordinates[count + 2] < epsilon &&
                   coordinates[count + 3] - expected_coordinates[count + 3] < epsilon &&
                   coordinates[count + 4] - expected_coordinates[count + 4] < epsilon );

            if( coordinates[count + 0] - expected_coordinates[count + 0] > epsilon ||
                coordinates[count + 1] - expected_coordinates[count + 1] > epsilon ||
                coordinates[count + 2] - expected_coordinates[count + 2] > epsilon ||
                coordinates[count + 3] - expected_coordinates[count + 3] > epsilon ||
                coordinates[count + 4] - expected_coordinates[count + 4] > epsilon  )
            {
                cout << endl << "MV_weights_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 5;
        }
    }

    cout << endl << "MV_weights_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
