// Author(s) : Dmitry Anisimov.
// We use an almost degenerate rectangle of length 1 and height 1.0e-200 and an exact data type
// in order to test Wachspress coordinates computed for some points from its interior.
// The test itself consists of generating some strictly interior points and then checking
// the constant and linear precision properties of obtained coordinates. Coordinates for the center
// point are also checked to be 1/4 exactly.

// Does not work with inexact kernel. We get division by zero due to W = 0.

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Wachspress_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::WP_coordinates_2<Polygon, Coordinate_vector> Wachspress_coordinates;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point vertices[4] = { Point(0, 0                                 ),
                                Point(1, 0                                 ),
                                Point(1, Scalar(1)/Scalar(pow(10.0, 200.0))),
                                Point(0, Scalar(1)/Scalar(pow(10.0, 200.0)))
                              };

    const Polygon rectangle(vertices, vertices + 4);

    Wachspress_coordinates wachspress_coordinates(rectangle);

    Coordinate_vector coordinates;

    const Point query_points[31] = { Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 300.0))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 280.0))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 260.0))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 240.0))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 220.0))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 210.0))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 200.8))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 200.5))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 200.2))) ,
                                     Point(Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 200.05))),

                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 300.0))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 280.0))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 260.0))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 240.0))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 220.0))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 210.0))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 200.8))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 200.5))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 200.2))) ,
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 200.05))),

                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 300.0))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 280.0))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 260.0))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 240.0))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 220.0))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 210.0))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 200.8))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 200.5))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 200.2))) ,
                                     Point(Scalar(1) - Scalar(1)/Scalar(pow(10.0, 300.0)), Scalar(1)/Scalar(pow(10.0, 200.05))),

                                     Point(Scalar(1)/Scalar(2),
                                           Scalar(1)/Scalar(2)*Scalar(rectangle.vertex(0).y() + rectangle.vertex(3).y()))
                                   };

    int count = 0;
    const Point zero(0, 0);

    for(int i = 0; i < 31; ++i) {
        const Output_type result = wachspress_coordinates.compute(query_points[i], coordinates);

        const Scalar coordinate_sum = coordinates[count + 0] +
                                      coordinates[count + 1] +
                                      coordinates[count + 2] +
                                      coordinates[count + 3] ;

        const Point linear_combination( rectangle.vertex(0).x()*coordinates[count + 0] +
                                        rectangle.vertex(1).x()*coordinates[count + 1] +
                                        rectangle.vertex(2).x()*coordinates[count + 2] +
                                        rectangle.vertex(3).x()*coordinates[count + 3] ,
                                        rectangle.vertex(0).y()*coordinates[count + 0] +
                                        rectangle.vertex(1).y()*coordinates[count + 1] +
                                        rectangle.vertex(2).y()*coordinates[count + 2] +
                                        rectangle.vertex(3).y()*coordinates[count + 3] );

        const Point difference(linear_combination.x() - query_points[i].x(), linear_combination.y() - query_points[i].y());

        if( (coordinate_sum != Scalar(1)) || (difference != zero) )
        {
            cout << endl << result.second << " WP_almost_degenerate_polygon_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 4;
    }

    const Scalar quater  = Scalar(1)/Scalar(4);
    if( coordinates[120] - quater != Scalar(0) ||
        coordinates[121] - quater != Scalar(0) ||
        coordinates[122] - quater != Scalar(0) ||
        coordinates[123] - quater != Scalar(0)  )
    {
        cout << endl << 1 << " WP_almost_degenerate_polygon_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "WP_almost_degenerate_polygon_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
