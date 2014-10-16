// Author(s) : Dmitry Anisimov.
// In this test we compute Mean Value coordinates for ~34400 strictly interior points with respect to
// a concave polygon with 10 vertices and check if they satisfy constant and linear precision properties.
// The chosen data type is inexact. The epsilon value is 1.0e-13.

// Works with an exact type, too.

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
    const Point vertices[10] = { Point(0                   , 0                   ),
                                 Point(Scalar(3) /Scalar(4), Scalar(1) /Scalar(2)),
                                 Point(Scalar(7) /Scalar(4), Scalar(-1)/Scalar(2)),
                                 Point(Scalar(11)/Scalar(4), Scalar(1) /Scalar(2)),
                                 Point(Scalar(7) /Scalar(2), 0                   ),
                                 Point(Scalar(7) /Scalar(2), 2                   ),
                                 Point(Scalar(11)/Scalar(4), Scalar(3) /Scalar(2)),
                                 Point(Scalar(7) /Scalar(4), Scalar(5) /Scalar(2)),
                                 Point(Scalar(3) /Scalar(4), Scalar(3) /Scalar(2)),
                                 Point(0                   , 2                   )
                               };

    const Polygon concave_polygon(vertices, vertices + 10);

    Mean_value_coordinates mean_value_coordinates(concave_polygon);

    Coordinate_vector coordinates;

    const Scalar step    = Scalar(1) / Scalar(100);
    const Scalar x_scale = Scalar(350);
    const Scalar y_scale = Scalar(100);

    int count = 0;
    const Scalar epsilon = Scalar(1) / Scalar(std::pow(10.0, 13.0));

    const Scalar limit_x = x_scale*step;
    const Scalar half = Scalar(1) / Scalar(2);
    const Scalar y_start = half + step;
    const Scalar limit_y = half + y_scale*step;

    for(Scalar x = step; x < limit_x; x += step) {
        for(Scalar y = y_start; y < limit_y; y += step) {
            const Point point(x, y);

            const Output_type result = mean_value_coordinates.compute(point, coordinates);

            const Scalar coordinate_sum = coordinates[count + 0] +
                                          coordinates[count + 1] +
                                          coordinates[count + 2] +
                                          coordinates[count + 3] +
                                          coordinates[count + 4] +
                                          coordinates[count + 5] +
                                          coordinates[count + 6] +
                                          coordinates[count + 7] +
                                          coordinates[count + 8] +
                                          coordinates[count + 9] ;

            const Point linear_combination( concave_polygon.vertex(0).x()*coordinates[count + 0] +
                                            concave_polygon.vertex(1).x()*coordinates[count + 1] +
                                            concave_polygon.vertex(2).x()*coordinates[count + 2] +
                                            concave_polygon.vertex(3).x()*coordinates[count + 3] +
                                            concave_polygon.vertex(4).x()*coordinates[count + 4] +
                                            concave_polygon.vertex(5).x()*coordinates[count + 5] +
                                            concave_polygon.vertex(6).x()*coordinates[count + 6] +
                                            concave_polygon.vertex(7).x()*coordinates[count + 7] +
                                            concave_polygon.vertex(8).x()*coordinates[count + 8] +
                                            concave_polygon.vertex(9).x()*coordinates[count + 9] ,
                                            concave_polygon.vertex(0).y()*coordinates[count + 0] +
                                            concave_polygon.vertex(1).y()*coordinates[count + 1] +
                                            concave_polygon.vertex(2).y()*coordinates[count + 2] +
                                            concave_polygon.vertex(3).y()*coordinates[count + 3] +
                                            concave_polygon.vertex(4).y()*coordinates[count + 4] +
                                            concave_polygon.vertex(5).y()*coordinates[count + 5] +
                                            concave_polygon.vertex(6).y()*coordinates[count + 6] +
                                            concave_polygon.vertex(7).y()*coordinates[count + 7] +
                                            concave_polygon.vertex(8).y()*coordinates[count + 8] +
                                            concave_polygon.vertex(9).y()*coordinates[count + 9] );

            const Point difference(linear_combination.x() - point.x(), linear_combination.y() - point.y());

            if( ((coordinate_sum - Scalar(1)) > epsilon) || difference.x() > epsilon || difference.y() > epsilon )
            {
                cout << endl << result.second << " MV_const_linear_precision_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 10;
        }
    }

    cout << endl << "MV_const_linear_precision_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
