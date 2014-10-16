// Author(s) : Dmitry Anisimov.
// In this test we compute Mean Value coordinates at some particular points,
// where the computation might break, with respect to a weakly convex polygon.
// We use inexact kernel, epsilon 1.0e-15, and we check if the resulting coordinate values are correct.

// Works with an exact kernel, too.

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mean_value_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

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
    const Point vertices[6] = { Point(0, 0                  ),
                                Point(1, 0                  ),
                                Point(1, Scalar(1)/Scalar(2)),
                                Point(1, 1                  ),
                                Point(0, 1                  ),
                                Point(0, Scalar(1)/Scalar(2))
                              };

    const Polygon weakly_convex_polygon(vertices, vertices + 6);

    Mean_value_coordinates mean_value_coordinates(weakly_convex_polygon);

    Coordinate_vector coordinates;

    const Point query_points[7] = { Point(Scalar(1) - (Scalar(1) / Scalar(std::pow(10.0, 300.0))),  Scalar(1) / Scalar(4)                     ),
                                    Point(Scalar(1) - (Scalar(1) / Scalar(std::pow(10.0, 300.0))),  Scalar(5) / Scalar(8)                     ),
                                    Point(Scalar(1) - (Scalar(1) / Scalar(std::pow(10.0, 300.0))),  Scalar(7) / Scalar(8)                     ),
                                    Point((Scalar(1) / Scalar(std::pow(10.0, 300.0)))            ,  Scalar(5) / Scalar(8)                     ),
                                    Point((Scalar(1) / Scalar(std::pow(10.0, 300.0)))            ,  Scalar(3) / Scalar(8)                     ),
                                    Point((Scalar(1) / Scalar(std::pow(10.0, 300.0)))            ,  Scalar(1) / Scalar(4)                     ),
                                    Point((Scalar(1) / Scalar(std::pow(10.0, 300.0)))            , (Scalar(1) / Scalar(std::pow(10.0, 300.0))))
                                   };

    const Scalar expected_coordinates[42] = { 0, Scalar(1)/Scalar(2), Scalar(1)/Scalar(2), 0, 0, 0,
                                              0, 0, Scalar(3)/Scalar(4), Scalar(1)/Scalar(4), 0, 0,
                                              0, 0, Scalar(1)/Scalar(4), Scalar(3)/Scalar(4), 0, 0,
                                              0, 0, 0, 0, Scalar(1)/Scalar(4), Scalar(3)/Scalar(4),
                                              Scalar(1)/Scalar(4), 0, 0, 0, 0, Scalar(3)/Scalar(4),
                                              Scalar(3)/Scalar(4), 0, 0, 0, 0, Scalar(1)/Scalar(4),
                                              1, 0, 0, 0, 0                                    , 0
                                            };

    int count = 0;
    const Scalar epsilon = Scalar(1) / Scalar(std::pow(10.0, 15.0));

    for(int i = 0; i < 5; ++i) {
        const Output_type result = mean_value_coordinates.compute(query_points[i], coordinates);

        if( CGAL::abs(coordinates[count + 0] - expected_coordinates[count + 0]) > epsilon ||
            CGAL::abs(coordinates[count + 1] - expected_coordinates[count + 1]) > epsilon ||
            CGAL::abs(coordinates[count + 2] - expected_coordinates[count + 2]) > epsilon ||
            CGAL::abs(coordinates[count + 3] - expected_coordinates[count + 3]) > epsilon ||
            CGAL::abs(coordinates[count + 4] - expected_coordinates[count + 4]) > epsilon ||
            CGAL::abs(coordinates[count + 5] - expected_coordinates[count + 5]) > epsilon  )
        {
            cout << endl << result.second << " MV_weakly_convex_polygon_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 6;
    }

    cout << endl << "MV_weakly_convex_polygon_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
