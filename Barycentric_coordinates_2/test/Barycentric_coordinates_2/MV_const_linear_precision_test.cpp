// Author: Dmitry Anisimov.
// In this test we compute mean value coordinates for ~34400 strictly interior points with respect to
// a concave polygon with 10 vertices and check if they satisfy constant and linear precision properties.
// The chosen data type is inexact. The epsilon value is 1.0e-13.

// Works with an exact type, too.

#include <cmath>
#include <cassert>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

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
    Point_vector vertices(10);

    vertices[0] = Point(0, 0);                                       vertices[1] = Point(Scalar(3) /Scalar(4), Scalar(1)/Scalar(2));
    vertices[2] = Point(Scalar(7) /Scalar(4), Scalar(-1)/Scalar(2)); vertices[3] = Point(Scalar(11)/Scalar(4), Scalar(1)/Scalar(2));
    vertices[4] = Point(Scalar(7) /Scalar(2), 0);                    vertices[5] = Point(Scalar(7) /Scalar(2), 2);
    vertices[6] = Point(Scalar(11)/Scalar(4), Scalar(3) /Scalar(2)); vertices[7] = Point(Scalar(7) /Scalar(4), Scalar(5)/Scalar(2));
    vertices[8] = Point(Scalar(3) /Scalar(4), Scalar(3) /Scalar(2)); vertices[9] = Point(0, 2);

    Mean_value_coordinates mean_value_coordinates(vertices.begin(), vertices.end());

    Coordinate_vector coordinates;

    const Scalar step    = Scalar(1) / Scalar(100);
    const Scalar x_scale = Scalar(350);
    const Scalar y_scale = Scalar(100);

    int count = 0;
    const Scalar epsilon = Scalar(1) / Scalar(std::pow(10.0, 8.0));

    const Scalar limit_x = x_scale*step;
    const Scalar half    = Scalar(1) / Scalar(2);
    const Scalar y_start = half + step;
    const Scalar limit_y = half + y_scale*step;

    for(Scalar x = step; x < limit_x; x += step) {
        for(Scalar y = y_start; y < limit_y; y += step) {
            const Point point(x, y);

            const Output_type result = mean_value_coordinates(point, coordinates);

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

            const Point linear_combination( vertices[0].x()*coordinates[count + 0] +
                                            vertices[1].x()*coordinates[count + 1] +
                                            vertices[2].x()*coordinates[count + 2] +
                                            vertices[3].x()*coordinates[count + 3] +
                                            vertices[4].x()*coordinates[count + 4] +
                                            vertices[5].x()*coordinates[count + 5] +
                                            vertices[6].x()*coordinates[count + 6] +
                                            vertices[7].x()*coordinates[count + 7] +
                                            vertices[8].x()*coordinates[count + 8] +
                                            vertices[9].x()*coordinates[count + 9] ,
                                            vertices[0].y()*coordinates[count + 0] +
                                            vertices[1].y()*coordinates[count + 1] +
                                            vertices[2].y()*coordinates[count + 2] +
                                            vertices[3].y()*coordinates[count + 3] +
                                            vertices[4].y()*coordinates[count + 4] +
                                            vertices[5].y()*coordinates[count + 5] +
                                            vertices[6].y()*coordinates[count + 6] +
                                            vertices[7].y()*coordinates[count + 7] +
                                            vertices[8].y()*coordinates[count + 8] +
                                            vertices[9].y()*coordinates[count + 9] );

            const Point difference(linear_combination.x() - point.x(), linear_combination.y() - point.y());

            assert( ((coordinate_sum - Scalar(1)) < epsilon) && difference.x() < epsilon && difference.y() < epsilon );

            if( ((coordinate_sum - Scalar(1)) > epsilon) || difference.x() > epsilon || difference.y() > epsilon )
            {
                cout << endl << "MV_const_linear_precision_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 10;
        }
    }

    cout << endl << "MV_const_linear_precision_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
