// Author(s) : Dmitry Anisimov.
// We use an almost degenerate triangle with one side of length 1 and an exact data type
// in order to test coordinates computed for some strictly interior points.
// The test itself consists of generating some strictly interior points and then checking
// the linear precision property of obtained coordinates. Some points close to the boundary
// up to 1.0e-300 are used, and coordinates for the center point are checked to be exactly 1/3.

// Does not work with inexact kernel! Get inconsistency when comparing difference with zero.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Triangle_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;
typedef Kernel::Triangle_2 Triangle;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Triangle, Vector_insert_iterator> Triangle_coordinates;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string; using std::pow;

int main()
{
    const Triangle triangle( Point(0, 0), Point(1, 0), Point(Scalar(-1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 200.0))) );

    Coordinate_vector coordinates;

    Triangle_coordinates triangle_coordinates(triangle);

    const Point query_points[12] = { Point(0, Scalar(1)/Scalar(pow(10.0, 201.0))),
                                     Point(0, Scalar(1)/Scalar(pow(10.0, 220.0))),
                                     Point(0, Scalar(1)/Scalar(pow(10.0, 240.0))),
                                     Point(0, Scalar(1)/Scalar(pow(10.0, 260.0))),
                                     Point(0, Scalar(1)/Scalar(pow(10.0, 280.0))),
                                     Point(0, Scalar(1)/Scalar(pow(10.0, 300.0))),
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 200.5))),
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 201.0))),
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 260.0))),
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 280.0))),
                                     Point(Scalar(1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 300.0))),
                                     Point( (Scalar(1)/Scalar(3))*Scalar(triangle.vertex(0).x() + triangle.vertex(1).x() + triangle.vertex(2).x()) ,
                                            (Scalar(1)/Scalar(3))*Scalar(triangle.vertex(0).y() + triangle.vertex(1).y() + triangle.vertex(2).y()) )
                                  };

    int count = 0;
    const Point zero(0, 0);
    for(int i = 0; i < 12; ++i) {
        const Output_type result = triangle_coordinates.compute(query_points[i], std::back_inserter(coordinates));

        const Point linear_combination( triangle.vertex(0).x()*coordinates[count + 0] +
                                        triangle.vertex(1).x()*coordinates[count + 1] +
                                        triangle.vertex(2).x()*coordinates[count + 2] ,
                                        triangle.vertex(0).y()*coordinates[count + 0] +
                                        triangle.vertex(1).y()*coordinates[count + 1] +
                                        triangle.vertex(2).y()*coordinates[count + 2] );

        const Point difference(linear_combination.x() - query_points[i].x(), linear_combination.y() - query_points[i].y());

        if(difference != zero)
        {
            cout << endl << result.second << " Almost_degenerate_triangle_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 3;
    }

    const Scalar third = Scalar(1)/Scalar(3);
    if( coordinates[33] - third != Scalar(0) ||
        coordinates[34] - third != Scalar(0) ||
        coordinates[35] - third != Scalar(0)  )
    {
        cout << endl << 1 << " Almost_degenerate_triangle_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "Almost_degenerate_triangle_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
