// Author(s) : Dmitry Anisimov.
// We use an almost degenerate segment of length 1.0e-200 and an exact data type
// in order to test coordinates computed for some points along its interior.
// The test itself consists of generating some strictly interior points and then checking
// the linear precision property of obtained coordinates. Coordinates for the center
// point are also checked to be 1/2 exactly.

// Does not work with inexact kernel! Get division by zero during the normalization step.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Segment_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT        Scalar;
typedef Kernel::Point_2   Point;
typedef Kernel::Segment_2 Segment;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Segment_coordinates_2<Segment, Vector_insert_iterator> Segment_coordinates;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string; using std::pow;

int main()
{
    const Segment segment( Point(0, 0), Point(Scalar(1)/Scalar(pow(10.0, 200.0)), 0) );

    Coordinate_vector coordinates;

    Segment_coordinates segment_coordinates(segment);

    const Point query_points[7] = { Point(Scalar(1)/Scalar(pow(10.0, 200.02)), 0),
                                    Point(Scalar(1)/Scalar(pow(10.0, 200.05)), 0),
                                    Point(Scalar(1)/Scalar(pow(10.0, 201.0)),  0),
                                    Point(Scalar(1)/Scalar(pow(10.0, 220.0)),  0),
                                    Point(Scalar(1)/Scalar(pow(10.0, 230.0)),  0),
                                    Point(Scalar(1)/Scalar(pow(10.0, 260.0)),  0),
                                    Point( (Scalar(1)/Scalar(2))*Scalar(segment.vertex(0).x() + segment.vertex(1).x()) ,
                                           (Scalar(1)/Scalar(2))*Scalar(segment.vertex(0).y() + segment.vertex(1).y()) )
                                  };

    int count = 0;
    const Point zero(0, 0);
    for(int i = 0; i < 7; ++i) {
        const Output_type result = segment_coordinates.compute(query_points[i], std::back_inserter(coordinates));

        const Point linear_combination( segment.vertex(0).x()*coordinates[count + 0] + segment.vertex(1).x()*coordinates[count + 1]  ,
                                        segment.vertex(0).y()*coordinates[count + 0] + segment.vertex(1).y()*coordinates[count + 1] );

        const Point difference(linear_combination.x() - query_points[i].x(), linear_combination.y() - query_points[i].y());

        if(difference != zero)
        {
            cout << endl << result.second << " Almost_degenerate_segment_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        count += 2;
    }

    const Scalar half = Scalar(1)/Scalar(2);
    if( coordinates[12] - half != Scalar(0) ||
        coordinates[13] - half != Scalar(0)  )
    {
        cout << endl << 1 << " Almost_degenerate_segment_test: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "Almost_degenerate_segment_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
