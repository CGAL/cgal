// Author: Dmitry Anisimov.
// We use an almost degenerate triangle with one side of length 1 and an exact data type
// in order to test coordinates computed for some strictly interior points.
// The test itself consists of generating some strictly interior points and then checking
// the linear precision property of obtained coordinates. Some points close to the boundary
// up to 1.0e-300 are used, and coordinates for the center point are checked to be exactly 1/3.

// Does not work with inexact kernel! Get inconsistency when comparing difference with zero.

#include <cmath>
#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string; using std::pow;

int main()
{
    const Point first_vertex  = Point(0, 0);
    const Point second_vertex = Point(1, 0);
    const Point third_vertex  = Point(Scalar(-1)/Scalar(2), Scalar(1)/Scalar(pow(10.0, 200.0)));

    Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);

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
                                     Point( (Scalar(1)/Scalar(3))*Scalar(first_vertex.x() + second_vertex.x() + third_vertex.x()) ,
                                            (Scalar(1)/Scalar(3))*Scalar(first_vertex.y() + second_vertex.y() + third_vertex.y()) )
                                  };

    Coordinate_vector coordinates;

    int count = 0;
    const Point zero(0, 0);
    for(int i = 0; i < 12; ++i) {
        const Output_type result = triangle_coordinates(query_points[i], std::back_inserter(coordinates));

        const Point linear_combination(  first_vertex.x()*coordinates[count + 0] +
                                        second_vertex.x()*coordinates[count + 1] +
                                         third_vertex.x()*coordinates[count + 2] ,
                                         first_vertex.y()*coordinates[count + 0] +
                                        second_vertex.y()*coordinates[count + 1] +
                                         third_vertex.y()*coordinates[count + 2] );

        const Point difference(linear_combination.x() - query_points[i].x(), linear_combination.y() - query_points[i].y());

        assert(difference == zero);

        if(difference != zero)
        {
            cout << endl << "Almost_degenerate_triangle_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 3;
    }

    const Scalar third = Scalar(1)/Scalar(3);

    assert(coordinates[33] - third == Scalar(0) && coordinates[34] - third == Scalar(0) && coordinates[35] - third == Scalar(0));

    if( coordinates[33] - third != Scalar(0) ||
        coordinates[34] - third != Scalar(0) ||
        coordinates[35] - third != Scalar(0)  )
    {
        cout << endl << "Almost_degenerate_triangle_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "Almost_degenerate_triangle_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
