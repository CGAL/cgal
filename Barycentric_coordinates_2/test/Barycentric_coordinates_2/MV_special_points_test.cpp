// Author(s) : Dmitry Anisimov.
// In this test we compute Mean Value coordinates at some particular points,
// where the computation might break. The used polygon is a concave polygon with 7 vertices.
// We also use inexact kernel and epsilon 1.0e-15.

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
    const Point vertices[7] = { Point(0                  , 0                  ),
                                Point(1                  , 1                  ),
                                Point(Scalar(7)/Scalar(4), Scalar(1)/Scalar(2)),
                                Point(Scalar(7)/Scalar(4), Scalar(5)/Scalar(2)),
                                Point(1                  , 2                  ),
                                Point(0                  , 3                  ),
                                Point(Scalar(1)/Scalar(2), Scalar(3)/Scalar(2))
                               };

    const Polygon concave_polygon(vertices, vertices + 7);

    Mean_value_coordinates mean_value_coordinates(concave_polygon);

    Coordinate_vector coordinates;

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

    int count = 0;
    const Scalar epsilon = Scalar(1) / Scalar(std::pow(10.0, 15.0));

    for(int i = 0; i < 11; ++i) {
        const Output_type result = mean_value_coordinates.compute(query_points[i], coordinates);

        assert(!std::isnan(CGAL::to_double(coordinates[count + 0])));
        assert(!std::isinf(CGAL::to_double(coordinates[count + 0])));

        assert(!std::isnan(CGAL::to_double(coordinates[count + 1])));
        assert(!std::isinf(CGAL::to_double(coordinates[count + 1])));

        assert(!std::isnan(CGAL::to_double(coordinates[count + 2])));
        assert(!std::isinf(CGAL::to_double(coordinates[count + 2])));

        assert(!std::isnan(CGAL::to_double(coordinates[count + 3])));
        assert(!std::isinf(CGAL::to_double(coordinates[count + 3])));

        assert(!std::isnan(CGAL::to_double(coordinates[count + 4])));
        assert(!std::isinf(CGAL::to_double(coordinates[count + 4])));

        assert(!std::isnan(CGAL::to_double(coordinates[count + 5])));
        assert(!std::isinf(CGAL::to_double(coordinates[count + 5])));

        assert(!std::isnan(CGAL::to_double(coordinates[count + 6])));
        assert(!std::isinf(CGAL::to_double(coordinates[count + 6])));

        const Scalar coordinate_sum = coordinates[count + 0] +
                                      coordinates[count + 1] +
                                      coordinates[count + 2] +
                                      coordinates[count + 3] +
                                      coordinates[count + 4] +
                                      coordinates[count + 5] +
                                      coordinates[count + 6] ;

        const Point linear_combination( concave_polygon.vertex(0).x()*coordinates[count + 0] +
                                        concave_polygon.vertex(1).x()*coordinates[count + 1] +
                                        concave_polygon.vertex(2).x()*coordinates[count + 2] +
                                        concave_polygon.vertex(3).x()*coordinates[count + 3] +
                                        concave_polygon.vertex(4).x()*coordinates[count + 4] +
                                        concave_polygon.vertex(5).x()*coordinates[count + 5] +
                                        concave_polygon.vertex(6).x()*coordinates[count + 6] ,
                                        concave_polygon.vertex(0).y()*coordinates[count + 0] +
                                        concave_polygon.vertex(1).y()*coordinates[count + 1] +
                                        concave_polygon.vertex(2).y()*coordinates[count + 2] +
                                        concave_polygon.vertex(3).y()*coordinates[count + 3] +
                                        concave_polygon.vertex(4).y()*coordinates[count + 4] +
                                        concave_polygon.vertex(5).y()*coordinates[count + 5] +
                                        concave_polygon.vertex(6).y()*coordinates[count + 6] );

        const Point difference(linear_combination.x() - query_points[i].x(), linear_combination.y() - query_points[i].y());

        if( ((coordinate_sum - Scalar(1)) > epsilon) || difference.x() > epsilon || difference.y() > epsilon )
        {
            cout << endl << result.second << " MV_special_points_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 7;
    }

    cout << endl << "MV_special_points_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
