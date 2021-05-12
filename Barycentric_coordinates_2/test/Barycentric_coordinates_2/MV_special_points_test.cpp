// Author: Dmitry Anisimov.
// In this test we compute mean value coordinates at some particular points,
// where the computation might break. The used polygon is a concave polygon with 7 vertices.
// We also use inexact kernel and epsilon = 1.0e-15.

// Works with an exact kernel, too.

#include <cmath>
#include <cassert>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>
#include <boost/math/special_functions/fpclassify.hpp>

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
    Point_vector vertices(7);

    vertices[0] = Point(0, 0);                                     vertices[1] = Point(1, 1);
    vertices[2] = Point(Scalar(7)/Scalar(4), Scalar(1)/Scalar(2)); vertices[3] = Point(Scalar(7)/Scalar(4), Scalar(5)/Scalar(2));
    vertices[4] = Point(1, 2);                                     vertices[5] = Point(0, 3);
    vertices[6] = Point(Scalar(1)/Scalar(2), Scalar(3)/Scalar(2));

    Mean_value_coordinates mean_value_coordinates(vertices.begin(), vertices.end());

    const Point query_points[11] = { Point(Scalar(1) + (Scalar(1) / Scalar(std::pow(10.0, 300.0)))            , Scalar(2) - (Scalar(1) / Scalar(std::pow(10.0, 300.0)))            ),
                                     Point(Scalar(1) + (Scalar(1) / Scalar(std::pow(10.0, 300.0)))            , Scalar(1) + (Scalar(1) / Scalar(std::pow(10.0, 300.0)))            ),
                                     Point(1                                                                  , Scalar(3) / Scalar(2)                                              ),
                                     Point(Scalar(5) / Scalar(4)                                              , Scalar(5) / Scalar(4)                                              ),
                                     Point(Scalar(5) / Scalar(4)                                              , Scalar(7) / Scalar(4)                                              ),
                                     Point(Scalar(3) / Scalar(2)                                              , Scalar(3) / Scalar(2)                                              ),

                                     Point(Scalar(7) / Scalar(4) - (Scalar(1) / Scalar(std::pow(10.0, 300.0))), Scalar(7) / Scalar(4) - (Scalar(1) / Scalar(std::pow(10.0, 300.0)))),
                                     Point(Scalar(7) / Scalar(4) - (Scalar(1) / Scalar(std::pow(10.0, 300.0))), Scalar(5) / Scalar(4) + (Scalar(1) / Scalar(std::pow(10.0, 300.0)))),

                                     Point(Scalar(3) / Scalar(4) - (Scalar(1) / Scalar(std::pow(10.0, 300.0))), Scalar(3) / Scalar(4) + (Scalar(1) / Scalar(std::pow(10.0, 300.0)))),
                                     Point(Scalar(3) / Scalar(4) - (Scalar(1) / Scalar(std::pow(10.0, 300.0))), Scalar(9) / Scalar(4) - (Scalar(1) / Scalar(std::pow(10.0, 300.0)))),
                                     Point(Scalar(1) / Scalar(2) + (Scalar(1) / Scalar(std::pow(10.0, 300.0))), Scalar(3) / Scalar(2)                                              )
                                   };

    Coordinate_vector coordinates;

    int count = 0;
    const Scalar epsilon = Scalar(1) / Scalar(std::pow(10.0, 15.0));

    for(int i = 0; i < 11; ++i) {
        const Output_type result = mean_value_coordinates(query_points[i], coordinates);

        assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 0])));
        assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 0])));

        assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 1])));
        assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 1])));

        assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 2])));
        assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 2])));

        assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 3])));
        assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 3])));

        assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 4])));
        assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 4])));

        assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 5])));
        assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 5])));

        assert(!boost::math::isnan(CGAL::to_double(coordinates[count + 6])));
        assert(!boost::math::isinf(CGAL::to_double(coordinates[count + 6])));

        const Scalar coordinate_sum = coordinates[count + 0] +
                                      coordinates[count + 1] +
                                      coordinates[count + 2] +
                                      coordinates[count + 3] +
                                      coordinates[count + 4] +
                                      coordinates[count + 5] +
                                      coordinates[count + 6] ;

        const Point linear_combination( vertices[0].x()*coordinates[count + 0] +
                                        vertices[1].x()*coordinates[count + 1] +
                                        vertices[2].x()*coordinates[count + 2] +
                                        vertices[3].x()*coordinates[count + 3] +
                                        vertices[4].x()*coordinates[count + 4] +
                                        vertices[5].x()*coordinates[count + 5] +
                                        vertices[6].x()*coordinates[count + 6] ,
                                        vertices[0].y()*coordinates[count + 0] +
                                        vertices[1].y()*coordinates[count + 1] +
                                        vertices[2].y()*coordinates[count + 2] +
                                        vertices[3].y()*coordinates[count + 3] +
                                        vertices[4].y()*coordinates[count + 4] +
                                        vertices[5].y()*coordinates[count + 5] +
                                        vertices[6].y()*coordinates[count + 6] );

        const Point difference(linear_combination.x() - query_points[i].x(), linear_combination.y() - query_points[i].y());

        assert( ((coordinate_sum - Scalar(1)) < epsilon) && difference.x() < epsilon && difference.y() < epsilon );

        if( ((coordinate_sum - Scalar(1)) > epsilon) || difference.x() > epsilon || difference.y() > epsilon )
        {
            cout << endl << "MV_special_points_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 7;
    }

    cout << endl << "MV_special_points_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
